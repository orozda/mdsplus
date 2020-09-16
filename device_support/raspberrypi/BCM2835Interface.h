#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <signal.h> 
#include <unistd.h>
#include <semaphore.h>

#define NUM_ADC_CHANNEL 4
#define NUM_DAC_CHANNEL 4

#define BCM2835_DAC_MAX_FREQ 10000.

#define INTERNAL 0
#define EXTERNAL 1


typedef struct adc_struct {
   short    *sampleData;
   uint64_t *sampleTime; 
   uint16_t  nSample;   
   uint64_t  periodMicoseconds;
   int      *stopAcq;
   sem_t    *mutex;
   uint8_t  dataValid;
   uint8_t  trigMode;
} t_adc_struct;

t_adc_struct bcm2835ADC;

extern "C" int spiConfig();
extern "C" int spiClose();
extern "C" int readData(short *sampleData, uint64_t *sampleTime, uint16_t numSample, uint64_t period);
extern "C" int bcn2835_readAndSaveAllChannels( int bufSize, int segmentSize, int numSamples, void *dataNidPtr, 
		        	               int clockNid, float timeIdx0, float period, void *treePtr, void *saveListPtr, void *stopAcq, 
                                               int shot, void *resampledNidPtr, int8_t trigMode);
extern "C" int getErrno();

extern "C" void getStopAcqFlag(void **stopAcq);
extern "C" void freeStopAcqFlag(void *stopAcq);
extern "C" void setStopAcqFlag(void *stopAcq);
extern "C" int bcm2835_write_ao(int ch, void *waveScaled, float gain, float offset, uint32_t numPoint);
extern "C" int bcm2835_start_ao(int fd, int trigMode);
extern "C" int bcm2835_stop_ao(int fd);

extern "C" void *bcm2835_generateWaves(void *arg);

typedef struct waves_struct {
   short *data;
   float offset;
   uint16_t numPoint;
} t_waves_struct;

typedef struct dac_struct {
   t_waves_struct waves[NUM_DAC_CHANNEL];
   int       stopGen;
   int       trigMode;
   sem_t     mutex;
   pthread_t thread;
} t_dac_struct;

t_dac_struct bcm2835DAC;

