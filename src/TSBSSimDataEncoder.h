#ifndef TSBSSIMDATAENCODER_H
#define TSBSSIMDATAENCODER_H

#include <vector>
#include <string>

#define SBS_MAX_ENCODER_WORDS 1024 // Max number of words the encoder can encode
#define SBS_NWORDS_MASK       0x3FFF
#define SBS_CHANNEL_MASK      0x1FF
#define SBS_TYPE_MASK         0x1FF
#define SBS_CHANNEL_FIRST_BIT 14
#define SBS_TYPE_FIRST_BIT    23

namespace SimEncoder {
  struct data {
    unsigned int channel;
  };

  struct adc_data : data {
    unsigned int integral;
  };

  struct fadc_data : adc_data {
    std::vector<unsigned int> samples;
  };

  struct tdc_data : data {
    std::vector<unsigned int> time;
    unsigned int getTime(unsigned int t) { return time[t]&0x7FFFFFFF; }
    unsigned int getEdge(unsigned int t) { return (time[t]>>31)&0x1; }
    // Note: bit 31 will be the edge (1 for trail, 0 for lead)
  };
};

class TSBSSimDataEncoder {
public:
  TSBSSimDataEncoder(const char *enc_name, unsigned short enc_id);
  virtual ~TSBSSimDataEncoder() {};

  // Encoders
  virtual bool EncodeADC(SimEncoder::adc_data data, unsigned int *enc_data,
      unsigned short &nwords) { return false; };
  virtual bool EncodeTDC(SimEncoder::tdc_data data, unsigned int *enc_data,
      unsigned short &nwords) { return false; };
  virtual bool EncodeFADC(SimEncoder::fadc_data data, unsigned int *enc_data,
      unsigned short &nwords) { return false; };
  // Decoders
  virtual bool DecodeADC(SimEncoder::adc_data &data,
      const unsigned int *enc_data,unsigned short nwords) { return false; }
  virtual bool DecodeTDC(SimEncoder::tdc_data &data,
      const unsigned int *enc_data,unsigned short nwords) { return false; };
  virtual bool DecodeFADC(SimEncoder::fadc_data &data,
      const unsigned int *enc_data,unsigned short nwords) { return false; }

  // Capabilities
  virtual bool IsADC() { return false; }
  virtual bool IsTDC() { return false; }
  virtual bool IsFADC() { return false; }

  unsigned short GetId() { return fEncId; }
  std::string GetName() { return fName; }

  static TSBSSimDataEncoder* GetEncoderByName(const char *enc_name);
  static TSBSSimDataEncoder* GetEncoder(unsigned short id);
  static unsigned int MakeBitMask(unsigned short bits);

  static unsigned int EncodeHeader(unsigned short type, unsigned short mult,
      unsigned int nwords);
  static void DecodeHeader(unsigned int hdr, unsigned short &type, unsigned short &ch,
      unsigned int &nwords);
  static unsigned short DecodeChannel(unsigned int hdr);
  static unsigned short DecodeType(unsigned int hdr);
  static unsigned short DecodeNwords(unsigned int hdr);

protected:
  std::string fName;
  unsigned short fEncId;

private:
  static std::vector<TSBSSimDataEncoder*> fEncoders;
};

// Generic TDC encoder
class TSBSSimTDCEncoder : public TSBSSimDataEncoder {
public:
  TSBSSimTDCEncoder(const char *enc_name, unsigned short enc_id,
      unsigned short bits, unsigned short edge_bit);
  virtual ~TSBSSimTDCEncoder() {};

  // Overloaded functions
  virtual bool EncodeTDC(SimEncoder::tdc_data data, unsigned int *enc_data,
      unsigned short &nwords);
  virtual bool DecodeTDC(SimEncoder::tdc_data &data,
      const unsigned int *enc_data,unsigned short nwords);
  virtual bool IsTDC() { return true; }

protected:
  unsigned short fBits;
  unsigned short fEdgeBit;
  unsigned short fBitMask;
};

// Generic ADC encoder
class TSBSSimADCEncoder : public TSBSSimDataEncoder {
public:
  TSBSSimADCEncoder(const char *enc_name, unsigned short enc_id,
      unsigned short bits);
  virtual ~TSBSSimADCEncoder() {};

  virtual bool EncodeADC(SimEncoder::adc_data data, unsigned int *enc_data,
      unsigned short &nwords);
  virtual bool DecodeADC(SimEncoder::adc_data &data,
      const unsigned int *enc_data,unsigned short nwords);
  virtual bool IsADC() { return true; }

protected:
  unsigned short fBits;
  unsigned short fBitMask;
};

// JLab FADC 250 in multi sample ADC mode
class TSBSSimFADC250Encoder : public TSBSSimADCEncoder {
public:
  TSBSSimFADC250Encoder(const char *enc_name, unsigned short enc_id);
  virtual ~TSBSSimFADC250Encoder() {};

  virtual bool EncodeFADC(SimEncoder::fadc_data data, unsigned int *enc_data,
      unsigned short &nwords);
  virtual bool DecodeFADC(SimEncoder::fadc_data &data,
      const unsigned int *enc_data,unsigned short nwords);
  virtual bool IsFADC() { return true; }

private:
  unsigned int EncodeSingleSample(unsigned int dat);
  void UnpackSamples(unsigned int enc_data,unsigned int *buff,
      bool *overflow, bool *valid);
};

#endif // TSBSSIMDATAENCODER_H
