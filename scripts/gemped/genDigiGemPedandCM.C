// This helper script will generate GEM pedestal and CM files needed for adding GEM ped and CM to digitized GEM ADC samples.

#include <string>
#include <map>
#include <iostream>
#include <fstream>
#include "TString.h"

using namespace std;

// Parameters to be read-in from the local database.
int db_nmodules {-1};
std::vector<int> db_modAPVmap {-1};

// Parameters to be read-in from the local MC database.
int db_mc_nmodules {-1};
std::vector<int> db_mc_modAPVmap {-1};

//module id, axis, position on axis
struct apvInfoGEM {
	int gemid, axis, pos, invert;

	void print() const {
	std::cout 
	<< "GEMId: " << gemid 
	<< " Axis: " << axis 
	<< " Pos: " << pos 
	<< " Invert: " << invert
	<< std::endl;
	};
	bool operator<(const apvInfoGEM& other) const {
	if (gemid != other.gemid) return gemid < other.gemid;
	if (axis != other.axis) return axis < other.axis;
	return pos < other.pos;
	}
};

//vtpcrate, fiber, adc_ch
struct apvInfoDAQ{
	//fill these values they id the APV
	int vtpcrate;
	int fiber; //NOTE: == mpd_id
	int adc_ch;
	// int invert;

	bool operator<(const apvInfoDAQ& other) const {
        if (vtpcrate != other.vtpcrate) return vtpcrate < other.vtpcrate;
        if (fiber != other.fiber) return fiber < other.fiber;
        return adc_ch < other.adc_ch;
    }

    void print() const {
        std::cout << "VTPcrate: " << vtpcrate 
                  << " Fiber: " << fiber 
                  << " ADC_ch: " << adc_ch 
                  << std::endl; 
	}
};

// Map with each APV's DAQ info (vtpcrate, fiber, adc_ch) being used as keys to get their values (gemid, axis, pos, invert) at the GEM module's end.
std::map<apvInfoDAQ, apvInfoGEM> db_apvInfoMap;
std::map<apvInfoGEM, apvInfoDAQ> db_mc_apvInfoMap; // For MC, the mapping is reversed: GEM info is used as keys to get DAQ info as values.

struct CommonMode {
  double mean{};
  double sigma{};
};

using cmmap = std::map<int, CommonMode>; // Map CM mean and sigma values to each APV number.
std::map<int, cmmap> trackercmmap; // Defining globally.

// 
struct Pedestal {
  double mean{};
  double rms{};
};

using pedmap = std::map<int, Pedestal>; // Map ped mean and rms value to each strip number (physical strip #, as used in the digitization).
std::map<int, pedmap> trackerpedmap; // Defining globally.



int readlocalDB ( const string& db_name );
int makeCMfile ( const string& cmfile_name, const string& cmfileout_name );
int makePedfile ( const string& pedfile_name, const string& pedfileout_name );

// Making the replay pedestal and common mode files for digitized GEM ADC samples.
int readlocalMCDB ( const string& mcdb_name );
int make_MCreplay_CMfile ( const string& mc_replay_cmfile_name );
int make_MCreplay_Pedfile ( const string& mc_replay_pedfile_name );


int genDigiGemPedandCM( const string outname_prepend, const string db_name, const string realdatcmfile_name, const string realdatpedfile_name, 
                        const bool makeMCreplayFiles = false, const string mcdb_name = "" ) {

	if ( readlocalDB(db_name) < 0 ) return -1; // Real real data local DB file and make APV and GEM module mapping.

    std::string digicmoutfile_name = outname_prepend + "_digiCM.dat";
    std::string digipedoutfile_name = outname_prepend + "_digiPed.dat";
    if ( makeCMfile(realdatcmfile_name,digicmoutfile_name) < 0 ) return -1;
    if ( makePedfile(realdatpedfile_name,digipedoutfile_name) < 0 ) return -1;

    if ( makeMCreplayFiles ) {
        if ( mcdb_name.empty() ) {
            std::cerr << "ERROR: MC DB file name is empty. Cannot generate MC replay files." << std::endl;
            return -1;
        }
        if ( readlocalMCDB(mcdb_name) < 0 ) return -1;
        std::string mc_replay_cmfile_name = outname_prepend + "_MCreplayDB_CM.dat";
        std::string mc_replay_pedfile_name = outname_prepend + "_MCreplayDB_Ped.dat";
        if ( make_MCreplay_CMfile(mc_replay_cmfile_name) < 0 ) return -1;
        if ( make_MCreplay_Pedfile(mc_replay_pedfile_name) < 0 ) return -1;
    }    

    return 0;
}



int readlocalDB(const std::string& db_name) {
    std::cout << "Reading Local Database file " << db_name << std::endl;
    std::ifstream m_dbFile(db_name);
    if (!m_dbFile.is_open()) {
        std::cerr << "ERROR: Could not open the DB file " << db_name << std::endl;
        return -1;
    }

    TString currentline;
    // First pass: find nmodules
    while (currentline.ReadLine(m_dbFile)) {
        if (!currentline.BeginsWith("#")) {
            TObjArray* tokens = currentline.Tokenize(" ");
            if (tokens->GetEntries() > 1) {
                TString skey = ((TObjString*)(*tokens)[0])->GetString();
                if (skey == "nmodules") {
                    db_nmodules = ((TObjString*)(*tokens)[1])->GetString().Atoi();
                }
            }
        }
    }

    if (db_nmodules <= 0) {
        std::cerr << "The number of modules not defined or set to 0 in the DB " << db_name << "!!!" << std::endl;
        return -1;
    }

    db_modAPVmap.clear();
    db_modAPVmap.resize(db_nmodules);

    m_dbFile.clear();
    m_dbFile.seekg(0); // Reset file pointer for second pass

    while (currentline.ReadLine(m_dbFile)) {
        if (!currentline.BeginsWith("#")) {
            TObjArray* tokens = currentline.Tokenize(" ");
            if (tokens->GetEntries() >= 1) {
                TString skey = ((TObjString*)(*tokens)[0])->GetString();
                for (int i = 0; i < db_nmodules; i++) {
                    if (skey == Form("m%i.chanmap", i)) {
                        while (currentline.ReadLine(m_dbFile)) {
                            if (currentline.BeginsWith("#")) break;
                            TObjArray* tokens_chanmap = currentline.Tokenize(" ");
                            if (tokens_chanmap->GetEntries() >= 9 &&
                                ((TObjString*)(*tokens_chanmap)[0])->GetString().IsDigit()) {
                                int vtpcrate = ((TObjString*)(*tokens_chanmap)[0])->GetString().Atoi();
                                int fiber    = ((TObjString*)(*tokens_chanmap)[2])->GetString().Atoi();
                                int adc_ch   = ((TObjString*)(*tokens_chanmap)[4])->GetString().Atoi();
                                int gemid    = i;
                                int pos      = ((TObjString*)(*tokens_chanmap)[6])->GetString().Atoi();
                                int invert   = ((TObjString*)(*tokens_chanmap)[7])->GetString().Atoi();
                                int axis     = ((TObjString*)(*tokens_chanmap)[8])->GetString().Atoi();

                                apvInfoGEM thisAPVgeminfo{gemid, axis, pos, invert};
                                apvInfoDAQ thisAPVdaqinfo{vtpcrate, fiber, adc_ch};
                                db_apvInfoMap[thisAPVdaqinfo] = thisAPVgeminfo;
                            } else break;
                            delete tokens_chanmap;
                        }
                    }                    
                }
            }
            if (tokens->GetEntries()>1){
               TString skey = ((TObjString*)(*tokens)[0])->GetString(); 
               for (int i=0; i<db_nmodules; i++){
                if(skey == Form("m%i.apvmap",i)){
                    int thismodAPVmap = ((TObjString*)(*tokens)[1])->GetString().Atoi();
                    db_modAPVmap[i] = thismodAPVmap;
                }
               }
            }
            delete tokens;
        }
    }    
    return 0;
}


int makeCMfile(const std::string& cmfile_name, const string& cmfileout_name) {
    std::cout << "Reading CM file " << cmfile_name << std::endl;
    std::ifstream m_cmfile(cmfile_name);
    if (!m_cmfile.is_open()) {
        std::cerr << "ERROR: could not open CM file " << cmfile_name << std::endl;
        return -1;
    }

    //std::map<int, cmmap> trackercmmap; // digiplane -> (pos -> CM info)
    trackercmmap.clear();

    TString currentline;

    while (currentline.ReadLine(m_cmfile)) {
        if (!currentline.BeginsWith("#")) {
            TObjArray* tokens = currentline.Tokenize(" ");
            if (tokens->GetEntries() >= 6 && ((TObjString*)(*tokens)[0])->GetString().IsDigit()) {
                int vtpcrate = ((TObjString*)(*tokens)[0])->GetString().Atoi();
                int fiber    = ((TObjString*)(*tokens)[2])->GetString().Atoi();
                int adc_ch   = ((TObjString*)(*tokens)[3])->GetString().Atoi();
                double cm_mean  = ((TObjString*)(*tokens)[4])->GetString().Atof();
                double cm_sigma = ((TObjString*)(*tokens)[5])->GetString().Atof();
                
                // Lookup GEM info
                apvInfoDAQ thisAPVdaqinfo{vtpcrate, fiber, adc_ch};
                auto it = db_apvInfoMap.find(thisAPVdaqinfo);
                if (it == db_apvInfoMap.end()) {
                    std::cerr << "Warning: APV not found in DB for crate=" << vtpcrate
                              << " fiber=" << fiber << " adc_ch=" << adc_ch << std::endl;
                    delete tokens;
                    continue;
                }

                int gemid = it->second.gemid;
                int axis  = it->second.axis;
                int pos   = it->second.pos;
                
                int digiplane_num = gemid * 2 + axis;
                trackercmmap[digiplane_num][pos] = CommonMode{cm_mean, cm_sigma};
            }
            delete tokens;
        }
    }

    std::ofstream m_cmdigifile(cmfileout_name);
    m_cmdigifile << "#For each GEM_plane#:\n";
    m_cmdigifile << "# APV# cm_mean cm_sigma\n" << '\n';

    for (const auto& [digiplane_num, apvcmmap] : trackercmmap) {
        m_cmdigifile << "GEM_plane# " << digiplane_num << "\n";
        for (const auto& [apvnum, cminfo] : apvcmmap) {
            m_cmdigifile << apvnum << " " << cminfo.mean << " " << cminfo.sigma << "\n";
        }
        m_cmdigifile << '\n';
    }

    m_cmdigifile.close();
    std::cout << "Output digi CM file " << cmfileout_name << " created.\n";
    return 0;
}


// APV 'raw strip' to 'physical strip' mapping.
namespace SBSGEM {
  enum GEMaxis_t { kUaxis=0, kVaxis };
  enum APVmap_t { kINFN=0, kUVA_XY, kUVA_UV, kMC };
}

constexpr int fN_APV25_CHAN {128};
std::array<std::vector<UInt_t>, 4 > APVMAP;

void InitAPVMAP(){
  APVMAP[SBSGEM::kINFN].resize(fN_APV25_CHAN);
  APVMAP[SBSGEM::kUVA_XY].resize(fN_APV25_CHAN);
  APVMAP[SBSGEM::kUVA_UV].resize(fN_APV25_CHAN);
  APVMAP[SBSGEM::kMC].resize(fN_APV25_CHAN);

  for( UInt_t i=0; i<fN_APV25_CHAN; i++ ){
    Int_t strip1 = 32*(i%4) + 8*(i/4) - 31*(i/16);
    Int_t strip2 = strip1 + 1 + strip1 % 4 - 5 * ( ( strip1/4 ) % 2 );
    Int_t strip3 = ( strip2 % 2 == 0 ) ? strip2/2 + 32 : ( (strip2<64) ? (63 - strip2)/2 : 127 + (65-strip2)/2 ); 
    APVMAP[SBSGEM::kINFN][i] = strip1; 
    APVMAP[SBSGEM::kUVA_XY][i] = strip2;
    APVMAP[SBSGEM::kUVA_UV][i] = strip3;
    APVMAP[SBSGEM::kMC][i] = i;
  } 
}

int GetStripNumber( UInt_t rawstrip, UInt_t pos, UInt_t invert, UInt_t modAPVmapping ){
  Int_t RstripNb = APVMAP[modAPVmapping][rawstrip];
  RstripNb = RstripNb + (127-2*RstripNb)*invert;
  Int_t RstripPos = RstripNb + 128*pos;

  // if( fIsMC ){
  //   return rawstrip + 128*pos;
  // }
  
  return RstripPos;
}

int makePedfile(const std::string& pedfile_name, const string& pedfileout_name){
    std::cout << "Reading pedestal file " << pedfile_name << std::endl;
    std::ifstream m_pedfile(pedfile_name);
    if (!m_pedfile.is_open()) {
        std::cerr << "ERROR: could not open ped file " << pedfile_name << std::endl;
        return -1;
    }

    InitAPVMAP();

    trackerpedmap.clear();
    TString currentline;
    int apv_vtpcrate = -1; 
    int apv_fiber    = -1;
    int apv_adc_ch   = -1;
    int apv_gemid    = -1;
    int apv_axis     = -1;
    int apv_pos      = -1;
    int apv_invert   = -1;
    int apv_digigemplane_num = -1;

    while(currentline.ReadLine(m_pedfile)) {
        if (!currentline.BeginsWith("#")) {
           TObjArray* tokens = currentline.Tokenize(" ");
           
           if (tokens->GetEntries() >= 5 && ((TObjString*)(*tokens)[0])->GetString()== "APV") {
            apv_vtpcrate = ((TObjString*)(*tokens)[1])->GetString().Atoi();
            apv_fiber = ((TObjString*)(*tokens)[3])->GetString().Atof();
            apv_adc_ch = ((TObjString*)(*tokens)[4])->GetString().Atof();

            // Lookup GEM info
            apvInfoDAQ thisAPVdaqinfo{apv_vtpcrate, apv_fiber, apv_adc_ch};
            auto it = db_apvInfoMap.find(thisAPVdaqinfo);
            if (it == db_apvInfoMap.end()) {
                std::cerr << "Warning: APV not found in DB for crate=" << apv_vtpcrate
                          << " fiber=" << apv_fiber << " adc_ch=" << apv_adc_ch << std::endl;
                delete tokens;
                continue;
            }

            apv_gemid  = it->second.gemid;
            apv_axis   = it->second.axis;
            apv_pos    = it->second.pos;
            apv_invert = it->second.invert;
            apv_digigemplane_num = apv_gemid*2 + apv_axis;
           }
           else if (tokens->GetEntries() >= 3 && ((TObjString*)(*tokens)[0])->GetString().IsDigit()) {
            int raw_stripnum =  ((TObjString*)(*tokens)[0])->GetString().Atoi();
            double pedmean = ((TObjString*)(*tokens)[1])->GetString().Atof();
            double pedrms = ((TObjString*)(*tokens)[2])->GetString().Atof();

            int physical_stripnum = GetStripNumber(raw_stripnum, apv_pos, apv_invert, db_modAPVmap.at(apv_gemid));
            trackerpedmap[apv_digigemplane_num][physical_stripnum] = Pedestal{pedmean,pedrms};
           }
           delete tokens;
        }
    }

    std::ofstream m_peddigifile(pedfileout_name);
    m_peddigifile << "#For each GEM_plane#:\n";
    m_peddigifile << "# Strip# ped_mean ped_rms\n" << '\n';

    for (const auto& [digiplane_num, apvpedmap] : trackerpedmap) {
        m_peddigifile << "GEM_plane# " << digiplane_num << "\n";
        for (const auto& [stripnum, pedinfo] : apvpedmap) {
            m_peddigifile << stripnum << " " << pedinfo.mean << " " << pedinfo.rms << "\n";
        }
        m_peddigifile << '\n';
    }

    m_peddigifile.close();
    std::cout << "Ouput digi ped file " << pedfileout_name << " created.\n";
    return 0;
}


int readlocalMCDB(const std::string& mcdb_name) {
    std::cout << "Reading Local MC Database file " << mcdb_name << std::endl;
    std::ifstream m_mcdbFile(mcdb_name);
    if (!m_mcdbFile.is_open()) {
        std::cerr << "ERROR: Could not open the MC DB file " << mcdb_name << std::endl;
        return -1;
    }

    TString currentline;
    // First pass: find nmodules
    while (currentline.ReadLine(m_mcdbFile)) {
        if (!currentline.BeginsWith("#")) {
            TObjArray* tokens = currentline.Tokenize(" ");
            if (tokens->GetEntries() > 1) {
                TString skey = ((TObjString*)(*tokens)[0])->GetString();
                if (skey == "nmodules") {
                    db_mc_nmodules = ((TObjString*)(*tokens)[1])->GetString().Atoi();
                }
            }
        }
    }

    if (db_mc_nmodules != db_nmodules) {
        std::cerr << "The number of modules defined in the MC DB " << mcdb_name << ", does not match the number defined in the real data DB!!!" << std::endl;
        return -1;
    }

    db_mc_modAPVmap.clear();
    db_mc_modAPVmap.resize(db_nmodules);

    m_mcdbFile.clear();
    m_mcdbFile.seekg(0); // Reset file pointer for second pass

    while (currentline.ReadLine(m_mcdbFile)) {
        if (!currentline.BeginsWith("#")) {
            TObjArray* tokens = currentline.Tokenize(" ");
            if (tokens->GetEntries() >= 1) {
                TString skey = ((TObjString*)(*tokens)[0])->GetString();
                for (int i = 0; i < db_nmodules; i++) {
                    if (skey == Form("m%i.chanmap", i)) {
                        while (currentline.ReadLine(m_mcdbFile)) {
                            if (currentline.BeginsWith("#")) break;
                            TObjArray* tokens_chanmap = currentline.Tokenize(" ");
                            if (tokens_chanmap->GetEntries() >= 9 &&
                                ((TObjString*)(*tokens_chanmap)[0])->GetString().IsDigit()) {
                                int vtpcrate = ((TObjString*)(*tokens_chanmap)[0])->GetString().Atoi();
                                int fiber    = ((TObjString*)(*tokens_chanmap)[2])->GetString().Atoi();
                                int adc_ch   = ((TObjString*)(*tokens_chanmap)[4])->GetString().Atoi();
                                int gemid    = i;
                                int pos      = ((TObjString*)(*tokens_chanmap)[6])->GetString().Atoi();
                                int invert   = ((TObjString*)(*tokens_chanmap)[7])->GetString().Atoi();
                                int axis     = ((TObjString*)(*tokens_chanmap)[8])->GetString().Atoi();

                                apvInfoGEM thisAPVgeminfo{gemid, axis, pos, invert};
                                apvInfoDAQ thisAPVdaqinfo{vtpcrate, fiber, adc_ch};
                                db_mc_apvInfoMap[thisAPVgeminfo] = thisAPVdaqinfo;
                            } else break;
                            delete tokens_chanmap;
                        }
                    }                    
                }
            }
            if (tokens->GetEntries()>1){
               TString skey = ((TObjString*)(*tokens)[0])->GetString(); 
               for (int i=0; i<db_nmodules; i++){
                if(skey == Form("m%i.apvmap",i)){
                    int thismodAPVmap = ((TObjString*)(*tokens)[1])->GetString().Atoi();
                    db_mc_modAPVmap[i] = thismodAPVmap;
                }
               }
            }
            delete tokens;
        }
    }    
    return 0;
}

int make_MCreplay_CMfile ( const string& mc_replay_cmfile_name ){

    std::ofstream mcreplay_cmfile(mc_replay_cmfile_name);
    if (  !mcreplay_cmfile.is_open() ) {
        std::cerr << "ERROR: could not open MC replay CM file " << mc_replay_cmfile_name << std::endl;
        return -1;
    }
    
    mcreplay_cmfile << "#format = crate, slot, mpd, adc_ch, CM mean, CM RMS:\n";

    for (const auto& [geminfo, daqinfo] : db_mc_apvInfoMap) {
        int gemid = geminfo.gemid;
        int axis  = geminfo.axis;
        int pos   = geminfo.pos;

        int digiplane_num = gemid * 2 + axis;

        auto cm_it = trackercmmap.find(digiplane_num);
        if (cm_it != trackercmmap.end()) {
            const auto& apvcmmap = cm_it->second;
            auto apv_cm_it = apvcmmap.find(pos);
            if (apv_cm_it != apvcmmap.end()) {
                const CommonMode& cminfo = apv_cm_it->second;
                mcreplay_cmfile << daqinfo.vtpcrate << " "
                                << daqinfo.fiber << " " // We use fiber ID for slot numer in the replay CM MC file.
                                << daqinfo.fiber << " "
                                << daqinfo.adc_ch << " "
                                << cminfo.mean << " "
                                << cminfo.sigma << "\n";
            } else {
                std::cerr << "Warning: No CM info found for GEM plane " << digiplane_num
                          << " and APV position " << pos << std::endl;
            }
        } else {
            std::cerr << "Warning: No CM info found for GEM plane " << digiplane_num << std::endl;
        }
    }

    mcreplay_cmfile.close();
    return 0;
}

int make_MCreplay_Pedfile ( const string& mc_replay_pedfile_name ){

    std::ofstream mcreplay_pedfile(mc_replay_pedfile_name);
    if (  !mcreplay_pedfile.is_open() ) {
        std::cerr << "ERROR: could not open MC replay ped file " << mc_replay_pedfile_name << std::endl;
        return -1;
    }
    
    mcreplay_pedfile << "# format = APV  crate  slotmpd_id  adc_ch followed by" << std::endl;
    mcreplay_pedfile << "# MC_APV_ channel_ number  pedestal mean  pedestal rms" << std::endl;

    for (const auto& [geminfo, daqinfo] : db_mc_apvInfoMap) {
        int gemid = geminfo.gemid;
        int axis  = geminfo.axis;
        int pos   = geminfo.pos;
        int firststrip_gemplane = pos * 128; // First physical strip number for this APV in the GEM plane.

        int digiplane_num = gemid * 2 + axis;

        // Print APV header line for each APV
        mcreplay_pedfile << "APV " << daqinfo.vtpcrate << " "
                         << daqinfo.fiber << " " // We use fiber ID for slot numer in the replay ped MC file.
                         << daqinfo.fiber << " "
                         << daqinfo.adc_ch << "\n";

        auto ped_it = trackerpedmap.find(digiplane_num);
        if (ped_it != trackerpedmap.end()) {            
            const auto& apvpedmap = ped_it->second;
            
            for ( int i = 0; i < 128; i++ ){
                int gemplane_stripnum = firststrip_gemplane + i; // Physical strip number in the GEM plane.
                auto strip_it = apvpedmap.find(gemplane_stripnum);
                if (strip_it != apvpedmap.end()) {
                    const Pedestal& pedinfo = strip_it->second; 
                    mcreplay_pedfile << i << " " // MC APV channel number == APV physical strip number.
                                     << pedinfo.mean << " "
                                     << pedinfo.rms << "\n";
                }
            }           
        } else {
            std::cerr << "Warning: No pedestal info found for GEM plane " << digiplane_num << std::endl;
        }
    }

    mcreplay_pedfile.close();
    return 0;
}