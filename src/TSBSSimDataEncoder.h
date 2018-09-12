#ifndef TSBSSIMDATAENCODER_H
#define TSBSSIMDATAENCODER_H

#include <vector>

#define SBS_MAX_ENCODER_WORDS 1024 // Max number of words the encoder can encode
#define SBS_NWORDS_MASK       0x3FFF
#define SBS_CHANNEL_MASK      0x1FF
#define SBS_TYPE_MASK         0x1FF
#define SBS_CHANNEL_FIRST_BIT 14
#define SBS_TYPE_FIRST_BIT    23

namespace SimEncoder {

  const unsigned short CaenV1190  = 0;
  const unsigned short FADC250    = 1;
  const unsigned short Lecroy1877 = 2;
  const unsigned short VETROC     = 3;
  const unsigned short F1TDC      = 4;

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
    unsigned int lead_time;
    unsigned int trail_time;
  };

  unsigned int EncodeHeader(unsigned short type, unsigned short mult,
      unsigned int nwords);
  void DecodeHeader(unsigned int hdr, unsigned short &type, unsigned short &ch,
      unsigned int &nwords);
  unsigned short DecodeChannel(unsigned int hdr);
  unsigned short DecodeType(unsigned int hdr);
  unsigned short DecodeNwords(unsigned int hdr);



  bool F1TDCEncode(tdc_data data, unsigned int *enc_data,
      unsigned short &nwords);
  bool F1TDCDecode(tdc_data &dat, const unsigned int *enc_data,
      unsigned short nwords);
  //bool F1TDCEncode(tdc_data data, std::vector<unsigned int> &enc_data);
  //unsigned int F1TDCDecode(tdc_data &data, unsigned int *enc_data );

  bool FADC250Encode(fadc_data data, unsigned int *enc_data,
      unsigned short &nwords);
  unsigned int FADC250EncodeSingleSample(unsigned int dat);
  bool FADC250Decode(fadc_data &dat, const unsigned int *enc_data,
      unsigned short nwords);
  void FADC250UnpackSamples(unsigned int enc_data,unsigned int *buff,
      bool *overflow, bool *valid);
};

#endif // TSBSSIMDATAENCODER_H
