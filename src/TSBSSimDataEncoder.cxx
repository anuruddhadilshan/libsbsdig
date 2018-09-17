#include "TSBSSimDataEncoder.h"
#include <TString.h>
#include <iostream>

// This one is static, so define it again here
std::vector<TSBSSimDataEncoder*> TSBSSimDataEncoder::fEncoders;

TSBSSimDataEncoder* TSBSSimDataEncoder::GetEncoderByName(
    const char *enc_name)
{
  if(fEncoders.empty()) { // First generate the list of known encoders!!
    unsigned short ids = 1;
    // TDCs
    fEncoders.push_back(new TSBSSimTDCEncoder("caen1190",ids++,19,26));
    fEncoders.push_back(new TSBSSimTDCEncoder("lecroy1877",ids++,16,16));
    fEncoders.push_back(new TSBSSimTDCEncoder("vetroc",ids++,16,26));
    fEncoders.push_back(new TSBSSimTDCEncoder("f1tdc",ids++,16,31));
    // ADCs
    fEncoders.push_back(new TSBSSimFADC250Encoder("fadc250",ids++));
    fEncoders.push_back(new TSBSSimADCEncoder("adc",ids++,12));
    fEncoders.push_back(new TSBSSimADCEncoder("lecroy1881",ids++,14));
    fEncoders.push_back(new TSBSSimADCEncoder("caen792",ids++,12));
  }

  TString name(enc_name);
  for(std::vector<TSBSSimDataEncoder*>::iterator it = fEncoders.begin();
      it != fEncoders.end(); it++) {
    if(name.CompareTo((*it)->GetName(),TString::kIgnoreCase)==0)
      return *it;
  }

  return 0;
}

TSBSSimDataEncoder* TSBSSimDataEncoder::GetEncoder(unsigned short id)
{
  for(std::vector<TSBSSimDataEncoder*>::iterator it = fEncoders.begin();
      it != fEncoders.end(); it++) {
    if((*it)->GetId() == id)
      return *it;
  }
  return 0;
}

TSBSSimDataEncoder::TSBSSimDataEncoder(const char *enc_name,
    unsigned short enc_id) : fName(enc_name), fEncId(enc_id)
{
}

TSBSSimTDCEncoder::TSBSSimTDCEncoder(const char *enc_name,
    unsigned short enc_id, unsigned short bits, unsigned short edge_bit)
  : TSBSSimDataEncoder(enc_name,enc_id), fBits(bits), fEdgeBit(edge_bit)
{
  fBitMask = MakeBitMask(fBits);
}

TSBSSimADCEncoder::TSBSSimADCEncoder(const char *enc_name,
    unsigned short enc_id, unsigned short bits)
  : TSBSSimDataEncoder(enc_name,enc_id), fBits(bits)
{
  fBitMask = MakeBitMask(fBits);
}

TSBSSimFADC250Encoder::TSBSSimFADC250Encoder(const char *enc_name,
    unsigned short enc_id) : TSBSSimADCEncoder(enc_name,enc_id,12)
{
}


unsigned int TSBSSimDataEncoder::MakeBitMask(unsigned short bits)
{
  unsigned int mask = 0;
  for(unsigned short b = 0; b < bits; b++) {
    mask |= 1<<b;
  }
  return mask;
}

bool TSBSSimADCEncoder::EncodeADC(SimEncoder::adc_data data,
    unsigned int *enc_data,unsigned short &nwords)
{
  nwords = 0;
  enc_data[nwords++] = data.integral&fBitMask;
  return (nwords>0);
}

bool TSBSSimTDCEncoder::EncodeTDC(SimEncoder::tdc_data data,
    unsigned int *enc_data,unsigned short &nwords)
{
  // Generic TDC encoder where the lowest n-bits are the time
  // measurement, and a single edge bit specifies either 0: lead, 1: trail
  // (We ignore the channel bit because it goes unused in this simulation)
  nwords = 0;
  for(std::vector<unsigned int>::iterator it = data.time.begin();
      it != data.time.end(); it++) {
    enc_data[nwords++] = ((*it)&fBitMask) | ((((*it)>>31)&0x1)<<fEdgeBit);
  }
  return (nwords>0);
}

bool TSBSSimFADC250Encoder::EncodeFADC(SimEncoder::fadc_data data,
    unsigned int *enc_data, unsigned short &nwords)
{
  // Word 1:
  //   (31) = 1 (<-- Skipping, not useful for Sim data)
  //   (30 – 27) = 4 (<-- I'm going to skip this, not sure it is useful)
  //   (26 – 23) = channel number (0 – 15) (<-- Also skipping)
  //   (22 – 12) = reserved (read as 0)
  //   (11 – 0) = window width (in number of samples)
  // Words 2 - N:
  //   (31) = 0
  //   (30) = reserved (read as 0)
  //   (29) = sample x not valid
  //   (28 – 16) = ADC sample x (includes overflow bit)
  //   (15 – 14) = reserved (read as 0)
  //   (13) = sample x + 1 not valid
  //   (12 – 0) = ADC sample x + 1 (includes overflow bit)
  nwords = 0;
  unsigned int nsamps = data.samples.size() >= 0xFFF ? 0xFFF :
    data.samples.size() ;
  enc_data[nwords++] = nsamps;
  unsigned int s = 0;
  unsigned int buff[2] = {0,0};
  for(s = 0; s < nsamps-1; s+=2) {
    buff[0] = EncodeSingleSample(data.samples[s]);
    buff[1] = EncodeSingleSample(data.samples[s+1]);
    enc_data[nwords++] = (buff[0]<<16) | buff[1];
  }
  if( s < nsamps ) { // Still have one more sample to process
    buff[0] = EncodeSingleSample(data.samples[s]);
    buff[1] = 0x2000; // Mark last sample in this two-sample word as not valid
    enc_data[nwords++] = (buff[0]<<16) | buff[1];
  }
  return (nwords>1);
}

bool TSBSSimADCEncoder::DecodeADC(SimEncoder::adc_data &data,
      const unsigned int *enc_data,unsigned short nwords)
{
  if(nwords>1)
    return false;
  unsigned short nread = 0;

  data.integral = enc_data[nread++]&fBitMask;
  return nread==nwords;
}

bool TSBSSimTDCEncoder::DecodeTDC(SimEncoder::tdc_data &data,
    const unsigned int *enc_data,unsigned short nwords)
{
  for(unsigned short n = 0; n < nwords; n++) {
    data.time.push_back(((enc_data[n]>>fEdgeBit)<<31) |
        (enc_data[n]&fBitMask));
  }
  return !data.time.empty();
}

bool TSBSSimFADC250Encoder::DecodeFADC(SimEncoder::fadc_data &data,
    const unsigned int *enc_data,unsigned short nwords)
{
  int nsamples = enc_data[0]&0xFFF;
  int nsamples_read = 0;

  unsigned int buff[2] = {0,0};
  bool overflow[2] = { false, false};
  bool valid[2] = {false, false};
  for(unsigned short n = 1; n < nwords; n++) {
    UnpackSamples(enc_data[n],buff,overflow,valid);
    for(short k = 0; k < 2; k++) {
      if(valid[k]) {
        data.samples.push_back(buff[k]);
        nsamples_read++;
      }
    }
  }

  if(nsamples_read != nsamples) {
    std::cerr << "Error, number of samples read (" << nsamples_read
      << "), does not match number of expected samples (" << nsamples
      << ")." << std::endl;
    return false;
  }
  return true;

}

unsigned int TSBSSimFADC250Encoder::EncodeSingleSample(unsigned int dat)
{
  if(dat&0xFFFFF000) { // Data too large, turn on overflow
    dat = 0x1FFF;
  }
  return dat&0x1FFF;
}

void TSBSSimFADC250Encoder::UnpackSamples(unsigned int enc_data,
    unsigned int *buff, bool *overflow, bool *valid)
{
  unsigned int tmp;
  for(int k = 0; k < 2; k++) {
    tmp = (k==0 ? (enc_data>>16) : enc_data)&0x3FFF;
    buff[k] = tmp&0xFFF;
    overflow[k] = tmp&0x1000;
    valid[k] = !(tmp&0x2000);
  }
}


unsigned int TSBSSimDataEncoder::EncodeHeader(unsigned short type,
    unsigned short mult, unsigned int nwords)
{
  // First word bits
  // 31-23: encoder type
  // 22-14: channel multiplier (to be converted to local channel by SimDecoder)
  // 13-0 : number of words that follow
  return ((type&SBS_TYPE_MASK)<<SBS_TYPE_FIRST_BIT) |
    ((mult&SBS_CHANNEL_MASK)<<SBS_CHANNEL_FIRST_BIT) |
    (nwords&SBS_NWORDS_MASK);
}

void TSBSSimDataEncoder::DecodeHeader(unsigned int hdr, unsigned short &type, unsigned short &ch,
    unsigned int &nwords)
{
  type = DecodeType(hdr);
  ch = DecodeChannel(hdr);
  nwords = DecodeNwords(hdr);
}

unsigned short TSBSSimDataEncoder::DecodeChannel(unsigned int hdr) {
  return (hdr>>SBS_CHANNEL_FIRST_BIT)&SBS_CHANNEL_MASK;
}

unsigned short TSBSSimDataEncoder::DecodeType(unsigned int hdr) {
  return (hdr>>SBS_TYPE_FIRST_BIT)&SBS_TYPE_MASK;
}

unsigned short TSBSSimDataEncoder::DecodeNwords(unsigned int hdr) {
  return hdr&SBS_NWORDS_MASK;
}

