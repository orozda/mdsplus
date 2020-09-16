#include "BCM2835Interface.h"
#include "AsyncStoreManager.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/ioctl.h>
#include <fcntl.h>
#include <time.h>
#include <getopt.h>
#include <bcm2835.h>
#include <sys/time.h>

#define ADC0 0
#define ADC1 1
#define DAC0 0
#define DAC1 1
#define CH0 0
#define CH1 1

// Rele' on pin GPIO17
//#define PIN_RELE RPI_V2_GPIO_P1_11
// Trigger on pin GPIO23
#define PIN_TRIG RPI_V2_GPIO_P1_16


#define PIN_RELE RPI_V2_GPIO_P1_11
#define rele_pin_init()  bcm2835_gpio_fsel(PIN_RELE, BCM2835_GPIO_FSEL_OUTP)
#define rele_pin_set()   bcm2835_gpio_write(PIN_RELE, HIGH)
#define rele_pin_clr()   bcm2835_gpio_write(PIN_RELE, LOW)

#define SIZEOF_SPI_BUFFER 3
char tx_buf[SIZEOF_SPI_BUFFER] = { (char)0x01, (char)0x80, (char)0x00 };
char rx_buf[SIZEOF_SPI_BUFFER];

int getErrno() { return errno; }

void mySleepMicroseconds(int64_t delay_micros)
{
   struct timeval now, pulse;
   int cycles, micros;
 
   while (1)
   {
      //cycles = 0;
      gettimeofday(&pulse, NULL);
      micros = 0;
      while (micros < delay_micros)
      {
         //++cycles;
         gettimeofday(&now, NULL);
         if (now.tv_sec > pulse.tv_sec) micros = 1000000L; else micros = 0;
         micros = micros + (now.tv_usec - pulse.tv_usec);

      }
      //printf("delay %d microseconds took %d cycles\n", delay_micros, cycles);
   }
}


/*inline*/ short  read_adc_raw(int chip, int channel) {
  bcm2835_spi_chipSelect(chip?BCM2835_SPI_CS1:BCM2835_SPI_CS0);
  tx_buf[1] = channel?0xC0:0x80;
  bcm2835_spi_transfernb(tx_buf, rx_buf, SIZEOF_SPI_BUFFER);
  return (rx_buf[1]<<8|rx_buf[2]);
}


int spiConfig()
{

  fprintf(stderr, "spiConfig Initialization\n");

  if (!bcm2835_init()) {
    fprintf(stderr, "bcm2835_init failed. Are you running as root??\n");
    return -1;
  }

  if (!bcm2835_spi_begin()) {
    fprintf(stderr, "bcm2835_spi_begin failedg. Are you running as root??\n");
    return -1;
  }

  // Set the RELE pin to be an output 
  rele_pin_init();
  
  // Set Trigger input line
  bcm2835_gpio_fsel(PIN_TRIG, BCM2835_GPIO_FSEL_INPT);
  bcm2835_gpio_set_pud(PIN_TRIG, BCM2835_GPIO_PUD_UP);

  bcm2835_spi_setBitOrder(BCM2835_SPI_BIT_ORDER_MSBFIRST);    // The default
  bcm2835_spi_setDataMode(BCM2835_SPI_MODE0);                 // The default
  //bcm2835_spi_setClockDivider(BCM2835_SPI_CLOCK_DIVIDER_256); // 0.9765625 MHz - 128 is OK, but 1.953125 MHz is out of specs
  bcm2835_spi_setClockDivider(BCM2835_SPI_CLOCK_DIVIDER_128); // 128 is OK, but 1.953125 MHz is out of specs
  bcm2835_spi_setChipSelectPolarity(BCM2835_SPI_CS0, LOW);    // the default
  bcm2835_spi_setChipSelectPolarity(BCM2835_SPI_CS1, LOW);    // the default
  bcm2835_spi_chipSelect(BCM2835_SPI_CS0);                    // The default

  return 1;
}

int spiClose()
{
    fprintf(stderr, "spiClose Initialization\n");
    rele_pin_clr();

    bcm2835_spi_end();
    bcm2835_close();
    return 1;
}

int waitForTrigger(int *stop)
{
   // Reads the current level on the specified pin and returns either HIGH or LOW (0 or 1).
     
     printf ("Wait for Trigger!!!\n") ;
     while( bcm2835_gpio_lev(PIN_TRIG) == 1 && !(*stop) )
	bcm2835_delayMicroseconds(100);
     if( !(*stop) )
     {
        printf ("TRIGGER!!!\n") ;
        return 1;
     }
     return 0;
}


void *readDataFun(void *arg)
{
    int64_t loopDt;
    int64_t deltaSleep;
    int64_t loopStart;
    struct timeval tv; 
    struct timespec req = {0};
    req.tv_sec = 0;

    t_adc_struct *st = (t_adc_struct *)arg;
    int16_t *sampleData[2];
    int idx=0;
   
    sampleData[0] = (int16_t*)calloc( st->nSample * NUM_ADC_CHANNEL, sizeof(int16_t));
    sampleData[1] = (int16_t*)calloc( st->nSample * NUM_ADC_CHANNEL, sizeof(int16_t));
   
     
    //Set rele
    rele_pin_set();

    //int extTrigger = false;
    //Wait for trigger
    if(st->trigMode)
        waitForTrigger(st->stopAcq);


    while( !*(st->stopAcq) )
    {
        if ( idx > 0 )
        {
             // pause 
             gettimeofday(&tv, NULL);
             loopDt = (((int64_t)tv.tv_sec * 1000000) + (int64_t)tv.tv_usec) - (int64_t)loopStart;
             deltaSleep = ((int64_t)st->periodMicoseconds-(int64_t)loopDt);              

           //fprintf(stderr, "periodMicoseconds[us]=%llu, loopDt[us]=%6llu, deltaSleep[us]=%6d\n", st->periodMicoseconds, loopDt, deltaSleep);

             if (deltaSleep > 0) {
                  bcm2835_delayMicroseconds(deltaSleep);
             } else
                  fprintf(stderr, "periodMicoseconds[us]=%llu, loopDt[us]=%6llu, deltaSleep[us]=%6d\n", st->periodMicoseconds, loopDt, deltaSleep);

        }

	for(int i = 0; i < st->nSample; i++)
	{
	    gettimeofday(&tv, NULL);
            loopStart=(((int64_t)tv.tv_sec * 1000000) + (int64_t)tv.tv_usec);
	    st->sampleTime[i]               = loopStart; 
	    sampleData[idx%2][i]               = read_adc_raw(ADC0,CH0);
	    sampleData[idx%2][i+st->nSample]   = read_adc_raw(ADC0,CH1);
	    sampleData[idx%2][i+2*st->nSample] = read_adc_raw(ADC1,CH0);
	    sampleData[idx%2][i+3*st->nSample] = read_adc_raw(ADC1,CH1);
  
            if ( i < st->nSample-1 )
            { 
                 // pause 
                 gettimeofday(&tv, NULL);
                 loopDt = (((int64_t)tv.tv_sec * 1000000) + (int64_t)tv.tv_usec) - (int64_t)loopStart;
                 deltaSleep = ((int64_t)st->periodMicoseconds-(int64_t)loopDt);              

                 //fprintf(stderr, "(%6d) periodMicoseconds[us]=%llu, loopDt[us]=%6llu, deltaSleep[us]=%6d\n", i, st->periodMicoseconds, loopDt, deltaSleep);

                 if (deltaSleep > 0) {
                      bcm2835_delayMicroseconds(deltaSleep);
                 }
                 continue;
              }
          }
          st->sampleData = sampleData[idx%2];
          //printf("UNLOCK interlock %d\n", *(st->stopAcq) );
          st->dataValid = 1;
          sem_post(st->mutex);
          idx++;
   }
   sleep(2);
   free(sampleData[0]);
   free(sampleData[1]);
   printf("Exit thread\n");  
   return NULL;
}

int readData(short *sampleData, uint64_t *sampleTime, uint16_t nSample, uint64_t period)
{
    int i;
    struct timeval tv; 
    int64_t startTime;

    struct timespec req = {0};
    req.tv_sec = 0;

    gettimeofday(&tv, NULL);
    startTime = (((int64_t)tv.tv_sec * 1000000) + (int64_t)tv.tv_usec);
    for(i = 0; i < nSample; i++)
    {
      gettimeofday(&tv, NULL);
      int64_t loopStart=(((int64_t)tv.tv_sec * 1000000) + (int64_t)tv.tv_usec);
      sampleTime[i]           = loopStart; 
      sampleData[i]           = read_adc_raw(ADC0,CH0);
      sampleData[i+nSample]   = read_adc_raw(ADC0,CH1);
      sampleData[i+2*nSample] = read_adc_raw(ADC1,CH0);
      sampleData[i+3*nSample] = read_adc_raw(ADC1,CH1);
      // pause 
      gettimeofday(&tv, NULL);
      int64_t loopDt = (((int64_t)tv.tv_sec * 1000000) + (int64_t)tv.tv_usec) - (int64_t)loopStart;
      int64_t deltaSleep = ((int64_t)period-(int64_t)loopDt);              

      //if(i%100 == 0)fprintf(stderr, "period[us]=%llu, loopDt[us]=%6llu, deltaSleep[us]=%6d\n", period, loopDt, deltaSleep);

      if (deltaSleep > 0) {
        bcm2835_delayMicroseconds(deltaSleep);
        //usleep(deltaSleep); minimum sleep 100us
      }

    }
    return nSample;
}

/***************************************/

int bcn2835_readAndSaveAllChannels( int bufSize, int segmentSize, int numSamples, void *dataNidPtr, 
			            int clockNid, float timeIdx0, float period, void *treePtr, void *saveListPtr, void *stopAcq, 
                                    int shot, void *resampledNidPtr, int8_t trigMode)

{

    int nChan = 4;
    int currSize = 0;
    int16_t *sampleData;
    uint64_t *sampleTime;
    double *time;
    uint16_t nSample;
    uint64_t periodMicoseconds;
    int currDataToRead = 0;
    bool transientRec 	= false;

    int slaveDataSamples;
    int __count = 0;
    int sampleToRead = 0;
    int readElem;

    SaveList *saveList = (SaveList *)saveListPtr; // Class to equeu data buffer to save in pulse file
    int *dataNid       = (int *)dataNidPtr;       // Channel node identifier
    int *resampledNid  = (int *)resampledNidPtr;  // Channel node identifier

    float dummyCalibCoeff[] = {1.,0.,0.,0.};

    char  *streamNames[nChan];
    float streamGains[nChan];
    float streamOffsets[nChan];
    //Dummy calibration coeffs
    float *coeffs[nChan];
    int numCoeffs[nChan];
    short *buffers_s[nChan]; 
    int readChanSmp[nChan];
    int bufReadChanSmp[nChan];
    float gains[nChan];
 

    printf("bcn2835_readAndSaveAllChannels \n");

    //Delete first all data nids
    for(int i = 0; i < nChan; i++)
    {
        try {
            TreeNode *currNode = new TreeNode(dataNid[i], (Tree *)treePtr);
            currNode->deleteData();
//Check if resampling
            try {
    	        Data *streamNameData = currNode->getExtendedAttribute("STREAM_NAME");
    	        streamNames[i] = streamNameData->getString();
    	        deleteData(streamNameData);
	        try {
		    Data *streamGainData = currNode->getExtendedAttribute("STREAM_GAIN");
		    streamGains[i] = streamGainData->getFloat();
	        }
	        catch(MdsException &exc) {streamGains[i] = 1;}
 	        try {
		    Data *streamOffsetData = currNode->getExtendedAttribute("STREAM_OFFSET");
		    streamOffsets[i] = streamOffsetData->getFloat();
	        }
	        catch(MdsException &exc) {streamOffsets[i] = 0;}
            }	
            catch(MdsException &exc) {streamNames[i] = NULL; streamGains[i] = 0; streamOffsets[i] = 0;}
            delete currNode;
	    if(resampledNid)
	    {
		  currNode = new TreeNode(resampledNid[i], (Tree *)treePtr);
		  currNode->deleteData();
		  delete currNode;
	    }
        }catch(MdsException &exc)
        {
            printf("Error deleting data nodes\n");
        }
    }

    for( int chan = 0; chan < nChan; chan++ )
    {
        buffers_s[chan] = new short[bufSize];
	bufReadChanSmp[chan] = bufSize;
	readChanSmp[chan] = 0;
        coeffs[chan] = dummyCalibCoeff;
        numCoeffs[chan] = 4;
        gains[chan] = 1.; 
    }
    
    nSample = bufSize;
    //sampleData = (int16_t*)calloc( nSample * NUM_ADC_CHANNEL, sizeof(int16_t));
    sampleTime = (uint64_t *)calloc(nSample, sizeof(uint64_t));
    periodMicoseconds = period * 1000000L;

    int streamFactor = (int)(0.1/period);
    if(bufSize > streamFactor)
	streamFactor = bufSize;
    else
    {
    	if(streamFactor % bufSize != 0)
	    streamFactor = (bufSize + 1)*(streamFactor / bufSize);
    }

    if( (*(int*)stopAcq) == 1)
        transientRec = true;
    else
        numSamples = segmentSize;

    (*(int*)stopAcq) = 0;

    int currSample = 0;


    pthread_t thread;
    //t_adc_struct bcm2835ADC;
    sem_t mutex;

    if (sem_init(&mutex, 0, 0) != 0)
    {
        printf("\n sem init failed\n");
        return 1;
    }

    bcm2835ADC.sampleData = sampleData;
    bcm2835ADC.sampleTime = sampleTime;
    bcm2835ADC.nSample    = bufSize;
    bcm2835ADC.periodMicoseconds = periodMicoseconds;
    bcm2835ADC.stopAcq    = (int *)stopAcq;
    bcm2835ADC.mutex      = &mutex;
    bcm2835ADC.trigMode   = trigMode;

    if( pthread_create(&thread, NULL, readDataFun, (void *)&bcm2835ADC) )
    {
        printf("Error creating acquisition thread  %d \n", (! *(int *)stopAcq));
        return 1;        
    }

    // cpu_set_t: This data set is a bitset where each bit represents a CPU.
    cpu_set_t cpuset;
    // CPU_ZERO: This macro initializes the CPU set set to be the empty set.
    CPU_ZERO(&cpuset);
    // CPU_SET: This macro adds cpu to the CPU set set.
    CPU_SET(0, &cpuset);

    if (pthread_setaffinity_np(thread, sizeof(cpu_set_t), &cpuset) != 0) 
    {
        printf("Error set thread affinity");
        return 1;        
    }

    while( (! *(uint8_t *)stopAcq ) == 1 )
    {   
            //Check if  has been acquired all sample 
            //and stop acquisition has not beeen asserted
            if( readChanSmp[0] == numSamples || (*(int*)stopAcq) )
            {
		if( transientRec  )
                {
                   (*(int*)stopAcq) = 1;
		   printf("Complete acquisition for all channels, read samples %d stop flag %d \n", readChanSmp[0], (*(int*)stopAcq) );            
                   break;
                }
                numSamples += numSamples; //continuous acquisition acquires a new segment
            }

            sampleToRead = numSamples - readChanSmp[0]; // Compute the sample to read for all channels
            if( sampleToRead < bufSize )
                currDataToRead = sampleToRead;
            else
                currDataToRead = bufSize;

            //currSample = readData(sampleData, sampleTime, currDataToRead, periodMicoseconds);
           
            //printf("WAIT for Buffer\n");
            bcm2835ADC.dataValid = 0;     
            sem_wait(&mutex);
            if ( bcm2835ADC.dataValid == 0 )
                break;
            //printf("UNLOCK for Buffer\n");
 
            sampleData = bcm2835ADC.sampleData;
            currSample = currDataToRead;  
            for( int ch = 0; ch <  NUM_ADC_CHANNEL; ch++ )
	    {
		 memcpy(buffers_s[ch], &sampleData[ch*currSample], currSample * sizeof(uint16_t));  
		 bufReadChanSmp[ch] = currSample;
	    }

            //Enqueue data to store in the pulse file
	    if( bufReadChanSmp[0] == currDataToRead )
	    {
  	         //Check if have been read more than required samples		
		 if ( readChanSmp[0] + bufReadChanSmp[0] > numSamples ) 
	 	        bufReadChanSmp[0] = numSamples - readChanSmp[0];

		 //Compute the number of samples to complete segment acquisition
	 	 sampleToRead = numSamples - readChanSmp[0];

		 for( int chan = 0; chan <  NUM_ADC_CHANNEL; chan++)
	 	 {
		      if(resampledNid)
		   	   saveList->addItem( buffers_s[chan], 
					      bufReadChanSmp[chan], sampleToRead, SHORT , segmentSize, 
					      readChanSmp[chan], dataNid[chan], clockNid, timeIdx0, treePtr, shot, streamFactor, streamNames[chan], 
					      streamGains[chan], streamOffsets[chan], period, gains[chan], coeffs[chan], numCoeffs[chan], resampledNid[chan]);
	              else
			   saveList->addItem( buffers_s[chan], 
					      bufReadChanSmp[chan], sampleToRead, SHORT, segmentSize, 
					      readChanSmp[chan], dataNid[chan], clockNid, timeIdx0, treePtr, shot, streamFactor, streamNames[chan], 
					      streamGains[chan], streamOffsets[chan], period, gains[chan], coeffs[chan], numCoeffs[chan]);

	              buffers_s[chan] = new short[bufSize];

		      //Update the number of samples read
		      readChanSmp[chan] += bufReadChanSmp[chan];
		      //Reset the the number of sample read for the next segment
		      bufReadChanSmp[chan] = 0;    
		 }
	   }
    } 

    printf("EXIT from bcn2835 ADC acquisition  %d \n", (! *(int *)stopAcq));

    pthread_join(thread,NULL);
    sem_destroy(&mutex);

  //free(sampleData);
    free(sampleTime);


    return 0;

}

void getStopAcqFlag(void **stopAcq)
{
    *stopAcq = (void *)malloc(sizeof(int));
}

void freeStopAcqFlag(void *stopAcq)
{
    free((char *)stopAcq);
}

void setStopAcqFlag(void *stopAcq)
{
    (*(int*)stopAcq) = 1;
    if(bcm2835ADC.mutex != 0)
       sem_post(bcm2835ADC.mutex);
}

/*********************************************************************************/
/*********************************************************************************/
/*********************************************************************************/
/****************************        DAC          ********************************/
/*********************************************************************************/
/*********************************************************************************/
/*********************************************************************************/

/*inline*/ void write_dac_raw(int chip, int channel, uint16_t raw) {
  uint16_t tx;
  uint16_t rx;
  bcm2835_spi_chipSelect(chip?BCM2835_SPI_CS1:BCM2835_SPI_CS0);
  if (channel) { raw = raw | 0xF000; } else { raw = (raw&0xFFF) | 0x3000; }
  // dac gain is 2
  raw = (raw &= ~(1 << 13));
  tx = (raw << 8) | (raw >> 8);
  bcm2835_spi_transfernb((char *)&tx, (char *)&rx, sizeof(tx));
  return;
}


int bcm2835_write_ao(int ch, void *waveScaled, float gain, float offset, uint32_t numPoint)
{
/*
   float offset = 0.;
   float gain = 4906. / 5.;
*/   
printf("bcm2835_write_ao chan %d\n",ch);

   bcm2835DAC.waves[ch].data = (short *) calloc( numPoint, sizeof(short) );
   for(int i = 0; i < numPoint; i++)
       bcm2835DAC.waves[ch].data[i] = (short)(gain * ((float*)waveScaled)[i] + offset );

   bcm2835DAC.waves[ch].numPoint = numPoint;
   bcm2835DAC.waves[ch].offset   = offset;
   bcm2835DAC.stopGen = 0;
   return numPoint;
}


int bcm2835_start_ao(int fd, int trigMode)
{
printf("bcm2835_start_ao mutex\n");
    if (sem_init(&bcm2835DAC.mutex , 0, 0) != 0)
    {
        printf("\n sem init failed\n");
        return 1;
    }

    bcm2835DAC.stopGen     = 0;
    bcm2835DAC.trigMode    = trigMode;

printf("bcm2835_start_ao thread\n");

    if( pthread_create(&bcm2835DAC.thread, NULL, bcm2835_generateWaves, (void *)&bcm2835DAC) )
    {
        printf("Error creating waveform generation thread  %d \n", bcm2835DAC.stopGen );
        return 1;        
    }

printf("bcm2835_start_ao affinity\n");

    // cpu_set_t: This data set is a bitset where each bit represents a CPU.
    cpu_set_t cpuset;
    // CPU_ZERO: This macro initializes the CPU set set to be the empty set.
    CPU_ZERO(&cpuset);
    // CPU_SET: This macro adds cpu to the CPU set set.
    CPU_SET(0, &cpuset);

    if (pthread_setaffinity_np(bcm2835DAC.thread, sizeof(cpu_set_t), &cpuset) != 0) 
    {
        printf("Error set thread affinity");
        return 1;        
    }

    return 0;
}

int bcm2835_stop_ao(int fd)
{
    bcm2835DAC.stopGen = 1;
    pthread_join(bcm2835DAC.thread, NULL);
printf("bcm2835_stop_ao\n");
    for( int ch=0; ch < NUM_DAC_CHANNEL; ch++)
    {
          if(bcm2835DAC.waves[ch].data != NULL)
              free(bcm2835DAC.waves[ch].data);
          bcm2835DAC.waves[ch].data=NULL;
          bcm2835DAC.waves[ch].numPoint=0;
    }
printf("bcm2835_stop_ao fine\n");
    return 0;
}

void *bcm2835_generateWaves(void *arg)
{
    int *outPtr;
    struct timeval tv; 
    t_dac_struct *wavesDAC = (t_dac_struct *)arg;
    uint64_t clock = 0;
    int64_t  period = 100;//10kHz waveform generation generation clock

    struct timespec req = {0};
    req.tv_sec = 0;

    //outPtr = (int *)malloc(sizeof(int));
    //*outPtr = 0;

    //Set rele
    rele_pin_set();

printf("bcm2835_generateWaves \n");

    int triggered;
    //Wait for trigger
    if( wavesDAC->trigMode == EXTERNAL )
        triggered = waitForTrigger(&wavesDAC->stopGen);
    
    //*outPtr = triggered;


    while( !wavesDAC->stopGen )
    {
        gettimeofday(&tv, NULL);
        int64_t loopStart=(((int64_t)tv.tv_sec * 1000000) + (int64_t)tv.tv_usec);


        //Generate waveform if trigger line is HIGH
       // Reads the current level on the specified pin and returns either HIGH or LOW (0 or 1).
       if ( bcm2835_gpio_lev(PIN_TRIG) == 0 )
       {
            write_dac_raw(DAC0,CH0,wavesDAC->waves[0].data[clock % wavesDAC->waves[0].numPoint]);
            write_dac_raw(DAC0,CH1,wavesDAC->waves[1].data[clock % wavesDAC->waves[1].numPoint]);
            write_dac_raw(DAC1,CH0,wavesDAC->waves[2].data[clock % wavesDAC->waves[2].numPoint]);
            write_dac_raw(DAC1,CH1,wavesDAC->waves[3].data[clock % wavesDAC->waves[3].numPoint]);

            if(clock%1000 == 0) printf("%llu %d %d %d %d\n", clock, wavesDAC->waves[0].data[clock % wavesDAC->waves[0].numPoint],   
                                         			wavesDAC->waves[1].data[clock % wavesDAC->waves[1].numPoint],
								wavesDAC->waves[2].data[clock % wavesDAC->waves[2].numPoint],
								wavesDAC->waves[3].data[clock % wavesDAC->waves[3].numPoint]);
        }
        else
        {
            write_dac_raw(DAC0,CH0,wavesDAC->waves[0].offset);
            write_dac_raw(DAC0,CH1,wavesDAC->waves[1].offset);
            write_dac_raw(DAC1,CH0,wavesDAC->waves[2].offset);
            write_dac_raw(DAC1,CH1,wavesDAC->waves[3].offset);
            if(clock%10000 == 0) printf("Trig Wave genration Disabled\n");
        }

        // pause 
        gettimeofday(&tv, NULL);
        int64_t loopDt = (((int64_t)tv.tv_sec * 1000000) + (int64_t)tv.tv_usec) - (int64_t)loopStart;
        int64_t deltaSleep = ((int64_t)period-(int64_t)loopDt);              

        //if(clock%1000 == 0)fprintf(stderr, "period[us]=%llu, loopDt[us]=%6llu, deltaSleep[us]=%6d\n", period, loopDt, deltaSleep);

        if (deltaSleep > 0) {
            //usleep(deltaSleep);
            bcm2835_delayMicroseconds(deltaSleep);
        }
        clock++;
    }
printf("fine bcm2835_generateWaves \n");
    //Set DAC to last wave point value
    write_dac_raw(DAC0,CH0,wavesDAC->waves[0].offset);
    write_dac_raw(DAC0,CH1,wavesDAC->waves[1].offset);
    write_dac_raw(DAC1,CH0,wavesDAC->waves[2].offset);
    write_dac_raw(DAC1,CH1,wavesDAC->waves[3].offset);
    //return (void *)outPtr;
    return NULL;
}

void setStopGenFlag(int state)
{
    bcm2835DAC.stopGen = state;
}

