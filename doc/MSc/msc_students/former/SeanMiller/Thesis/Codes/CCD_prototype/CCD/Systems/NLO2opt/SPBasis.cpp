#include "SPBasis.hpp"
#include <algorithm>
#include <cstdlib>
#include <vector>
#include <omp.h>

double SPBasis::referenceEnergy(){
  double Eref = 0.;
  // for(std::size_t i = 0; i < this->nParticles-1; i++){
  //   for(std::size_t j = i+1; j < this->nParticles; j++){
  for(std::size_t i = 0; i < this->nParticles; i++){
    for(std::size_t j = 0; j < this->nParticles; j++){
      Eref += 0.5*this->calcTBME(i,j,i,j);
    } // end j
  } // end i
  for(std::size_t i = 0; i < this->nParticles; i++){
    Eref += this->spEnergy[i];
  }
  return Eref;
}

void SPBasis::rotateSpEnergiesToNormalOrdered(){
  double ei,ea;
  for(std::size_t i = 0; i < this->nParticles; i++){
    ei = this->spEnergy[i];
    for(std::size_t j = 0; j < this->nParticles; j++){
      ei += this->calcTBME(i,j,i,j);
    }
    this->spEnergy[i] = ei;
  }
  for(std::size_t a = this->nParticles; a < this->nSpstates; a++){
    ea = this->spEnergy[a];
    for(std::size_t j = 0; j < this->nParticles; j++){
      ea += this->calcTBME(a,j,a,j);
    }
    this->spEnergy[a] = ea;
  }
}


void SPBasis::setUpChannelDims(){
  // first count how many unique channels there are
  // could potentially nuke these loops and just overallocate
  // upper bound is known

  this->chanDims = new ChannelDims[this->nChannels];
  this->chanModDims = new ChannelDims[this->nChannels];
  this->threeBodyChanDims = new ThreeBodyChannelDims[this->nSpstates];
  // delete this at some point
  this->channelIndices = new size_t[this->nChannels];
  this->inverseChannelIndices = new size_t[this->nChannels];
  this->modChannelIndices = new size_t[this->nChannels];
  this->inverseModChannelIndices = new size_t[this->nChannels];
  this->spChannelIndices = new size_t[this->nSpstates];
  this->inverseSPChannelIndices = new size_t[this->nSpstates];
  #pragma omp parallel
  {
    // initialize these to the identity
    #pragma omp for
    for(std::size_t i = 0; i < this->nChannels; i++){
      this->channelIndices[i] = i;
      this->inverseChannelIndices[i] = i;
      this->modChannelIndices[i] = i;
      this->inverseModChannelIndices[i] = i;
    }

    #pragma omp for
    for(std::size_t i = 0; i < this->nSpstates; i++){
      this->spChannelIndices[i] = i;
      this->inverseSPChannelIndices[i] = i;
    }

    #pragma omp for
    for( std::size_t ichan = 0; ichan < this->nChannels; ichan++){
      this->chanDims[ichan].ppDim = 0;
      this->chanDims[ichan].hhDim = 0;
      this->chanDims[ichan].hpDim = 0;
      this->chanDims[ichan].phDim = 0;
    }

    // set up dimension
    #pragma omp for collapse(2)
    for(std::size_t p = 0; p < this->nSpstates; p++){
      for(std::size_t q = 0; q < this->nSpstates; q++){
        std::size_t iTB = this->TBchanIndexFunction(p,q);

        if( p < this->nParticles ) {
          if( q < this->nParticles ) {
            #pragma omp atomic
            this->chanDims[iTB].hhDim++;
          } else {
            #pragma omp atomic
            this->chanDims[iTB].hpDim++;
          }
        } else {
          if( q < this->nParticles ){
            #pragma omp atomic
            this->chanDims[iTB].phDim++;
          } else {
            #pragma omp atomic
            this->chanDims[iTB].ppDim++;
          }
        }
      }
    }

    #pragma omp for
    for( std::size_t ichan = 0; ichan < this->nChannels; ichan++){
      this->chanModDims[ichan].ppDim = 0;
      this->chanModDims[ichan].hhDim = 0;
      this->chanModDims[ichan].hpDim = 0;
      this->chanModDims[ichan].phDim = 0;
    }


    // set up dimension
    #pragma omp for collapse(2)
    for(std::size_t p = 0; p < this->nSpstates; p++){
      for(std::size_t q = 0; q < this->nSpstates; q++){
        std::size_t iTB = this->TBmodChanIndexFunction(p,q);
        if( p < this->nParticles ) {
          if( q < this->nParticles ) {
            #pragma omp atomic
            this->chanModDims[iTB].hhDim++;
          } else {
            #pragma omp atomic
            this->chanModDims[iTB].hpDim++;
          }
        } else {
          if( q < this->nParticles ){
            #pragma omp atomic
            this->chanModDims[iTB].phDim++;
          } else {
            #pragma omp atomic
            this->chanModDims[iTB].ppDim++;
          }
        }
      }
    }

    #pragma omp for
    for(std::size_t i = 0; i < this->nParticles; i++){
      // find Dim
      this->threeBodyChanDims[i].hppDim = 0;
      this->threeBodyChanDims[i].phhDim = 0;
    }

    #pragma omp for collapse(3)
    for(std::size_t j = 0; j < this->nParticles; j++){
      for(std::size_t a = this->nParticles; a < this->nSpstates; a++){
        for(std::size_t b = this->nParticles; b < this->nSpstates; b++){
          std::size_t sp_i;
          if(spIndexExists_from3Body(j,a,b) && (sp_i = spIndex_from3Body(j,a,b)) < this->nParticles ){
            #pragma omp atomic
            this->threeBodyChanDims[sp_i].hppDim++;
          }
        }
      }
    }

    // find Dim
    #pragma omp for
    for(std::size_t a = this->nParticles; a < this->nSpstates; a++){
      this->threeBodyChanDims[a].hppDim = 0;
      this->threeBodyChanDims[a].phhDim = 0;
    } // end a loop

    #pragma omp for collapse(3)
    for(std::size_t b = this->nParticles; b < this->nSpstates; b++){
      for(std::size_t i = 0; i < this->nParticles; i++){
        for(std::size_t j = 0; j < this->nParticles; j++){
          std::size_t sp_a;
          if(spIndexExists_from3Body(b,i,j) && (sp_a = spIndex_from3Body(b,i,j)) >= this->nParticles ){
            #pragma omp atomic
            this->threeBodyChanDims[sp_a].phhDim++;
          }
        } // end j
      } // end i
    } // end b
  }

  this->maxChannelSize_h_hpp = 0;
  for(std::size_t i = 0; i < this->nParticles; i++){
    this->maxChannelSize_h_hpp = std::max(this->maxChannelSize_h_hpp, this->threeBodyChanDims[i].hppDim);
  }

  this->maxChannelSize_phh_p = 0;
  for(std::size_t a = this->nParticles; a < this->nSpstates; a++){
    this->maxChannelSize_phh_p = std::max(this->maxChannelSize_phh_p, this->threeBodyChanDims[a].phhDim);
  }

  this->maxChannelSize_hhhh = 0;
  this->maxChannelSize_hhhp = 0;
  this->maxChannelSize_hhpp = 0;
  this->maxChannelSize_hpph = 0;
  this->maxChannelSize_hphp = 0;
  this->maxChannelSize_hppp = 0;
  this->maxChannelSize_pppp = 0;

  for( std::size_t ichan = 0; ichan < this->nChannels; ichan++){
    std::size_t ppDim = this->chanModDims[ichan].ppDim;
    std::size_t hhDim = this->chanModDims[ichan].hhDim;
    std::size_t hpDim = this->chanModDims[ichan].hpDim;
    std::size_t phDim = this->chanModDims[ichan].phDim;

    this->maxChannelSize_hhhh = std::max(this->maxChannelSize_hhhh, hhDim * hhDim);
    this->maxChannelSize_hhhp = std::max(this->maxChannelSize_hhhp, hhDim * hpDim);
    this->maxChannelSize_hhpp = std::max(this->maxChannelSize_hhpp, hhDim * ppDim);
    this->maxChannelSize_hpph = std::max(this->maxChannelSize_hpph, hpDim * phDim);
    this->maxChannelSize_hphp = std::max(this->maxChannelSize_hphp, hpDim * hpDim);
    this->maxChannelSize_hppp = std::max(this->maxChannelSize_hppp, hpDim * ppDim);
    this->maxChannelSize_pppp = std::max(this->maxChannelSize_pppp, ppDim * ppDim);
  }
} // end setUpChannelDims


struct ChannelMapMaker{
  std::vector<PQBundle> ppMap;
  std::vector<PQBundle> hhMap;
  std::vector<PQBundle> hpMap;
  std::vector<PQBundle> phMap;
};

struct ThreeBodyChannelMapMaker{
  std::vector<PQRBundle> phhMap;
  std::vector<PQRBundle> hppMap;
};

void SPBasis::setUpChannelMaps(){
  this->chanModMaps = new ChannelMaps[this->nChannels];
  this->chanMaps = new ChannelMaps[this->nChannels];
  this->threeBodyChanMaps = new ThreeBodyChannelMaps[this->nSpstates];

  #pragma omp parallel
  {
    // allocate maps
    #pragma omp for schedule(guided)
    for(std::size_t ichan=0; ichan < this->nChannels; ichan++){
      this->chanMaps[ichan].ppMap = new PQBundle[this->chanDims[ichan].ppDim];
      this->chanMaps[ichan].hhMap = new PQBundle[this->chanDims[ichan].hhDim];
      this->chanMaps[ichan].hpMap = new PQBundle[this->chanDims[ichan].hpDim];
      this->chanMaps[ichan].phMap = new PQBundle[this->chanDims[ichan].phDim];

      this->chanDims[ichan].ppDim = 0;
      this->chanDims[ichan].hhDim = 0;
      this->chanDims[ichan].hpDim = 0;
      this->chanDims[ichan].phDim = 0;
    }

    // set up maps
    #pragma omp for collapse(2) schedule(guided)
    for(std::size_t p = 0; p < this->nSpstates; p++){
      for(std::size_t q = 0; q < this->nSpstates; q++){
        std::size_t iTB = this->TBchanIndexFunction(p,q);
        std::size_t atomicIdx;

        if( p < this->nParticles ) {
          if( q < this->nParticles ) {
            #pragma omp atomic capture
            atomicIdx = this->chanDims[iTB].hhDim++;
            this->chanMaps[iTB].hhMap[atomicIdx].p = p;
            this->chanMaps[iTB].hhMap[atomicIdx].q = q;
          } else {
            #pragma omp atomic capture
            atomicIdx = this->chanDims[iTB].hpDim++;
            this->chanMaps[iTB].hpMap[atomicIdx].p = p;
            this->chanMaps[iTB].hpMap[atomicIdx].q = q;
          }
        } else {
          if( q < this->nParticles ){
            #pragma omp atomic capture
            atomicIdx = this->chanDims[iTB].phDim++;
            this->chanMaps[iTB].phMap[atomicIdx].p = p;
            this->chanMaps[iTB].phMap[atomicIdx].q = q;
          } else {
            #pragma omp atomic capture
            atomicIdx = this->chanDims[iTB].ppDim++;
            this->chanMaps[iTB].ppMap[atomicIdx].p = p;
            this->chanMaps[iTB].ppMap[atomicIdx].q = q;
          }
        }

      }
    }

    #pragma omp for schedule(dynamic)
    for(std::size_t ichan=0; ichan < this->nChannels; ichan++){
      std::sort(this->chanMaps[ichan].ppMap, this->chanMaps[ichan].ppMap + this->chanDims[ichan].ppDim);
      std::sort(this->chanMaps[ichan].phMap, this->chanMaps[ichan].phMap + this->chanDims[ichan].phDim);
      std::sort(this->chanMaps[ichan].hpMap, this->chanMaps[ichan].hpMap + this->chanDims[ichan].hpDim);
      std::sort(this->chanMaps[ichan].hhMap, this->chanMaps[ichan].hhMap + this->chanDims[ichan].hhDim);
    }
    //////////////////////////////////////////////////////////////////////////////////////////////
    // Modified Channels now!!!!!!!
    /////////////////////////////////////////////////////////////////////////////////////////////


    // allocate maps
    #pragma omp for schedule(guided)
    for(std::size_t ichan=0; ichan < this->nChannels; ichan++){
      this->chanModMaps[ichan].ppMap = new PQBundle[this->chanModDims[ichan].ppDim];
      this->chanModMaps[ichan].hhMap = new PQBundle[this->chanModDims[ichan].hhDim];
      this->chanModMaps[ichan].hpMap = new PQBundle[this->chanModDims[ichan].hpDim];
      this->chanModMaps[ichan].phMap = new PQBundle[this->chanModDims[ichan].phDim];

      this->chanModDims[ichan].ppDim = 0;
      this->chanModDims[ichan].hhDim = 0;
      this->chanModDims[ichan].hpDim = 0;
      this->chanModDims[ichan].phDim = 0;

    }


    // set up maps
    #pragma omp for collapse(2) schedule(guided)
    for(std::size_t p = 0; p < this->nSpstates; p++){
      for(std::size_t q = 0; q < this->nSpstates; q++){
        std::size_t iTB = this->TBmodChanIndexFunction(p,q);
        std::size_t atomicIdx;

        if( p < this->nParticles ) {
          if( q < this->nParticles ) {
            #pragma omp atomic capture
            atomicIdx = this->chanModDims[iTB].hhDim++;
            this->chanModMaps[iTB].hhMap[atomicIdx].p = p;
            this->chanModMaps[iTB].hhMap[atomicIdx].q = q;
          } else {
            #pragma omp atomic capture
            atomicIdx = this->chanModDims[iTB].hpDim++;
            this->chanModMaps[iTB].hpMap[atomicIdx].p = p;
            this->chanModMaps[iTB].hpMap[atomicIdx].q = q;
          }
        } else {
          if( q < this->nParticles ){
            #pragma omp atomic capture
            atomicIdx = this->chanModDims[iTB].phDim++;
            this->chanModMaps[iTB].phMap[atomicIdx].p = p;
            this->chanModMaps[iTB].phMap[atomicIdx].q = q;
          } else {
            #pragma omp atomic capture
            atomicIdx = this->chanModDims[iTB].ppDim++;
            this->chanModMaps[iTB].ppMap[atomicIdx].p = p;
            this->chanModMaps[iTB].ppMap[atomicIdx].q = q;
          }
        }
      }
    }

    #pragma omp for schedule(dynamic)
    for(std::size_t ichan=0; ichan < this->nChannels; ichan++){
      std::sort(this->chanModMaps[ichan].ppMap, this->chanModMaps[ichan].ppMap + this->chanModDims[ichan].ppDim);
      std::sort(this->chanModMaps[ichan].phMap, this->chanModMaps[ichan].phMap + this->chanModDims[ichan].phDim);
      std::sort(this->chanModMaps[ichan].hpMap, this->chanModMaps[ichan].hpMap + this->chanModDims[ichan].hpDim);
      std::sort(this->chanModMaps[ichan].hhMap, this->chanModMaps[ichan].hhMap + this->chanModDims[ichan].hhDim);
    }

    /////////////////////////////////////////
    // 3-body channels
    /////////////////////////////////////////

    #pragma omp for schedule(guided)
    for(std::size_t i = 0; i < this->nParticles; i++){
      // allocate
      this->threeBodyChanMaps[i].hppMap = new PQRBundle[this->threeBodyChanDims[i].hppDim];
      // Now use dim as an index in the next section
      this->threeBodyChanDims[i].hppDim = 0;
    }

    // set
    #pragma omp for collapse(3) schedule(guided)
    for(std::size_t j = 0; j < this->nParticles; j++){
      for(std::size_t a = this->nParticles; a < this->nSpstates; a++){
        for(std::size_t b = this->nParticles; b < this->nSpstates; b++){
          std::size_t sp_i;
          std::size_t atomicIdx;
          if(spIndexExists_from3Body(j,a,b) &&
             (sp_i = spIndex_from3Body(j,a,b)) < this->nParticles ){
            #pragma omp atomic capture
            atomicIdx = this->threeBodyChanDims[sp_i].hppDim++;
            this->threeBodyChanMaps[sp_i].hppMap[atomicIdx].p = j;
            this->threeBodyChanMaps[sp_i].hppMap[atomicIdx].q = a;
            this->threeBodyChanMaps[sp_i].hppMap[atomicIdx].r = b;
          }
        } // end b loop
      } // end a loop
    } // end j loop



    // allocate
    #pragma omp for schedule(guided)
    for(std::size_t a = this->nParticles; a < this->nSpstates; a++){
      this->threeBodyChanMaps[a].phhMap = new PQRBundle[this->threeBodyChanDims[a].phhDim];
      // now use Dim as an index in next loop
      this->threeBodyChanDims[a].phhDim = 0;
    }
    // set

    #pragma omp for collapse(3) schedule(guided)
    for(std::size_t b = this->nParticles; b < this->nSpstates; b++){
      for(std::size_t i = 0; i < this->nParticles; i++){
        for(std::size_t j = 0; j < this->nParticles; j++){
          std::size_t sp_a;
          std::size_t atomicIdx;
          if(spIndexExists_from3Body(b,i,j) &&
             (sp_a = spIndex_from3Body(b,i,j)) >= this->nParticles ){
            #pragma omp atomic capture
            atomicIdx = this->threeBodyChanDims[sp_a].phhDim++;
            this->threeBodyChanMaps[sp_a].phhMap[atomicIdx].p = b;
            this->threeBodyChanMaps[sp_a].phhMap[atomicIdx].q = i;
            this->threeBodyChanMaps[sp_a].phhMap[atomicIdx].r = j;
          }
        }
      }
    }

    #pragma omp for schedule(dynamic)
    for(std::size_t ichan=0; ichan < this->nSpstates; ichan++){
      if(ichan < this->nParticles){
        std::sort(this->threeBodyChanMaps[ichan].hppMap, this->threeBodyChanMaps[ichan].hppMap + this->threeBodyChanDims[ichan].hppDim);
      }else{
        std::sort(this->threeBodyChanMaps[ichan].phhMap, this->threeBodyChanMaps[ichan].phhMap + this->threeBodyChanDims[ichan].phhDim);
      }
    }
  }
  /*
  ChannelMapMaker **chanMapMakers;
  ChannelMapMaker **chanModMapMakers;
  ThreeBodyChannelMapMaker **threeBodyChanMapMakers;

  #pragma omp parallel
  {
    int thread_num = omp_get_thread_num();
    int num_threads = omp_get_num_threads();
    #pragma omp single
    {
      chanMapMakers = new ChannelMapMaker*[num_threads];
      chanModMapMakers = new ChannelMapMaker*[num_threads];
      threeBodyChanMapMakers = new ThreeBodyChannelMapMaker*[num_threads];
    }
    chanMapMakers[thread_num] = new ChannelMapMaker[this->nChannels];
    chanModMapMakers[thread_num] = new ChannelMapMaker[this->nChannels];
    threeBodyChanMapMakers[thread_num] = new ThreeBodyChannelMapMaker[this->nSpstates];

    // set up maps
    #pragma omp for collapse(2) schedule(static, 1)
    for(std::size_t p = 0; p < this->nSpstates; p++){
      for(std::size_t q = 0; q < this->nSpstates; q++){
        std::size_t iTB = this->TBchanIndexFunction(p,q);

        if( p < this->nParticles ) {
          if( q < this->nParticles ) {
            chanMapMakers[thread_num][iTB].hhMap.push_back({p, q});
          } else {
            chanMapMakers[thread_num][iTB].hpMap.push_back({p, q});
          }
        } else {
          if( q < this->nParticles ){
            chanMapMakers[thread_num][iTB].phMap.push_back({p, q});
          } else {
            chanMapMakers[thread_num][iTB].ppMap.push_back({p, q});
          }
        }
      }
    }


    // allocate maps
    #pragma omp for schedule(dynamic)
    for(std::size_t ichan=0; ichan < this->nChannels; ichan++){
      this->chanMaps[ichan].hhMap = new PQBundle[this->chanDims[ichan].hhDim];
      this->chanMaps[ichan].hpMap = new PQBundle[this->chanDims[ichan].hpDim];
      this->chanMaps[ichan].phMap = new PQBundle[this->chanDims[ichan].phDim];
      this->chanMaps[ichan].ppMap = new PQBundle[this->chanDims[ichan].ppDim];

      std::size_t ihh = 0;
      std::size_t ihp = 0;
      std::size_t iph = 0;
      std::size_t ipp = 0;
      for(int t = 0; t < num_threads; t++){
        for(std::vector<PQBundle>::iterator pq = chanMapMakers[t][ichan].hhMap.begin(); pq != chanMapMakers[t][ichan].hhMap.end(); ++pq) {
          this->chanMaps[ichan].hhMap[ihh] = *pq;
          ihh++;
        }
        for(std::vector<PQBundle>::iterator pq = chanMapMakers[t][ichan].hpMap.begin(); pq != chanMapMakers[t][ichan].hpMap.end(); ++pq) {
          this->chanMaps[ichan].hpMap[ihp] = *pq;
          ihp++;
        }
        for(std::vector<PQBundle>::iterator pq = chanMapMakers[t][ichan].phMap.begin(); pq != chanMapMakers[t][ichan].phMap.end(); ++pq) {
          this->chanMaps[ichan].phMap[iph] = *pq;
          iph++;
        }
        for(std::vector<PQBundle>::iterator pq = chanMapMakers[t][ichan].ppMap.begin(); pq != chanMapMakers[t][ichan].ppMap.end(); ++pq) {
          this->chanMaps[ichan].ppMap[ipp] = *pq;
          ipp++;
        }
      }
    }

    //////////////////////////////////////////////////////////////////////////////////////////////
    // Modified Channels now!!!!!!!
    /////////////////////////////////////////////////////////////////////////////////////////////


    #pragma omp for collapse(2) schedule(static, 1)
    for(std::size_t p = 0; p < this->nSpstates; p++){
      for(std::size_t q = 0; q < this->nSpstates; q++){
        std::size_t iTB = this->TBmodChanIndexFunction(p,q);

        if( p < this->nParticles ) {
          if( q < this->nParticles ) {
            chanModMapMakers[thread_num][iTB].hhMap.push_back({p, q});
          } else {
            chanModMapMakers[thread_num][iTB].hpMap.push_back({p, q});
          }
        } else {
          if( q < this->nParticles ){
            chanModMapMakers[thread_num][iTB].phMap.push_back({p, q});
          } else {
            chanModMapMakers[thread_num][iTB].ppMap.push_back({p, q});
          }
        }
      }
    }

    // allocate maps
    #pragma omp for schedule(dynamic)
    for(std::size_t ichan=0; ichan < this->nChannels; ichan++){
      this->chanModMaps[ichan].hhMap = new PQBundle[this->chanModDims[ichan].hhDim];
      this->chanModMaps[ichan].hpMap = new PQBundle[this->chanModDims[ichan].hpDim];
      this->chanModMaps[ichan].phMap = new PQBundle[this->chanModDims[ichan].phDim];
      this->chanModMaps[ichan].ppMap = new PQBundle[this->chanModDims[ichan].ppDim];

      std::size_t ihh = 0;
      std::size_t ihp = 0;
      std::size_t iph = 0;
      std::size_t ipp = 0;
      for(int t = 0; t < num_threads; t++){
        for(std::vector<PQBundle>::iterator pq = chanModMapMakers[t][ichan].hhMap.begin(); pq != chanModMapMakers[t][ichan].hhMap.end(); ++pq) {
          this->chanModMaps[ichan].hhMap[ihh] = *pq;
          ihh++;
        }
        for(std::vector<PQBundle>::iterator pq = chanModMapMakers[t][ichan].hpMap.begin(); pq != chanModMapMakers[t][ichan].hpMap.end(); ++pq) {
          this->chanModMaps[ichan].hpMap[ihp] = *pq;
          ihp++;
        }
        for(std::vector<PQBundle>::iterator pq = chanModMapMakers[t][ichan].phMap.begin(); pq != chanModMapMakers[t][ichan].phMap.end(); ++pq) {
          this->chanModMaps[ichan].phMap[iph] = *pq;
          iph++;
        }
        for(std::vector<PQBundle>::iterator pq = chanModMapMakers[t][ichan].ppMap.begin(); pq != chanModMapMakers[t][ichan].ppMap.end(); ++pq) {
          this->chanModMaps[ichan].ppMap[ipp] = *pq;
          ipp++;
        }
      }
    }


    /////////////////////////////////////////
    // 3-body channels
    /////////////////////////////////////////


    // set
    #pragma omp for collapse(3) schedule(static, 1)
    for(std::size_t j = 0; j < this->nParticles; j++){
      for(std::size_t a = this->nParticles; a < this->nSpstates; a++){
        for(std::size_t b = this->nParticles; b < this->nSpstates; b++){
          std::size_t sp_i;
          std::size_t atomicIdx;
          if(spIndexExists_from3Body(j,a,b) &&
             (sp_i = spIndex_from3Body(j,a,b)) < this->nParticles ){
            threeBodyChanMapMakers[thread_num][sp_i].hppMap.push_back({j, a, b});
          }
        }
      }
    }

    #pragma omp for collapse(3) schedule(static, 1)
    for(std::size_t b = this->nParticles; b < this->nSpstates; b++){
      for(std::size_t i = 0; i < this->nParticles; i++){
        for(std::size_t j = 0; j < this->nParticles; j++){
          std::size_t sp_a;
          std::size_t atomicIdx;
          if(spIndexExists_from3Body(b,i,j) &&
             (sp_a = spIndex_from3Body(b,i,j)) >= this->nParticles ){
            threeBodyChanMapMakers[thread_num][sp_a].phhMap.push_back({b, i, j});
          }
        }
      }
    }

    #pragma omp for schedule(dynamic)
    for(std::size_t ichan=0; ichan < this->nSpstates; ichan++){
      if(ichan < this->nParticles){
        this->threeBodyChanMaps[ichan].hppMap = new PQRBundle[this->threeBodyChanDims[ichan].hppDim];

        std::size_t ihpp = 0;

        for(int t = 0; t < num_threads; t++){
          for(std::vector<PQRBundle>::iterator pqr = threeBodyChanMapMakers[t][ichan].hppMap.begin(); pqr != threeBodyChanMapMakers[t][ichan].hppMap.end(); ++pqr) {
            this->threeBodyChanMaps[ichan].hppMap[ihpp] = *pqr;
            ihpp++;
          }
        }
      } else {
        this->threeBodyChanMaps[ichan].phhMap = new PQRBundle[this->threeBodyChanDims[ichan].phhDim];

        std::size_t iphh = 0;

        for(int t = 0; t < num_threads; t++){
          for(std::vector<PQRBundle>::iterator pqr = threeBodyChanMapMakers[t][ichan].phhMap.begin(); pqr != threeBodyChanMapMakers[t][ichan].phhMap.end(); ++pqr) {
            this->threeBodyChanMaps[ichan].phhMap[iphh] = *pqr;
            iphh++;
          }
        }
      }
    }
    delete[] chanMapMakers[thread_num];
    delete[] chanModMapMakers[thread_num];
    delete[] threeBodyChanMapMakers[thread_num];
  }

  delete[] chanMapMakers;
  delete[] chanModMapMakers;
  delete[] threeBodyChanMapMakers;
  */
} // end setUpChannelMaps


void SPBasis::setUpInverseChannelIndicesAndPermuteDims(){
  ChannelDims * chanDimsTemp = new ChannelDims[this->nChannels];
  ChannelDims * chanModDimsTemp = new ChannelDims[this->nChannels];
  ThreeBodyChannelDims * threeBodyChanDimsTemp = new ThreeBodyChannelDims[this->nSpstates];

  for(std::size_t i = 0; i < this->nChannels; i++){
    // set up inverse map
    inverseChannelIndices[ channelIndices[i] ] = i;
    inverseModChannelIndices[ modChannelIndices[i] ] = i;

    // store dims in temps
    chanDimsTemp[i].ppDim = this->chanDims[i].ppDim;
    chanDimsTemp[i].hhDim = this->chanDims[i].hhDim;
    chanDimsTemp[i].phDim = this->chanDims[i].phDim;
    chanDimsTemp[i].hpDim = this->chanDims[i].hpDim;

    chanModDimsTemp[i].ppDim = this->chanModDims[i].ppDim;
    chanModDimsTemp[i].hhDim = this->chanModDims[i].hhDim;
    chanModDimsTemp[i].phDim = this->chanModDims[i].phDim;
    chanModDimsTemp[i].hpDim = this->chanModDims[i].hpDim;
  }

  for(std::size_t i = 0; i < this->nSpstates; i++){
    inverseSPChannelIndices[ spChannelIndices[i] ] = i;

    threeBodyChanDimsTemp[i].phhDim = this->threeBodyChanDims[i].phhDim;
    threeBodyChanDimsTemp[i].hppDim = this->threeBodyChanDims[i].hppDim;
  }

  for(std::size_t i = 0; i < this->nChannels; i++){
    this->chanDims[i].ppDim = chanDimsTemp[ this->channelIndices[i] ].ppDim;
    this->chanDims[i].hhDim = chanDimsTemp[ this->channelIndices[i] ].hhDim;
    this->chanDims[i].phDim = chanDimsTemp[ this->channelIndices[i] ].phDim;
    this->chanDims[i].hpDim = chanDimsTemp[ this->channelIndices[i] ].hpDim;

    this->chanModDims[i].ppDim = chanModDimsTemp[ this->modChannelIndices[i] ].ppDim;
    this->chanModDims[i].hhDim = chanModDimsTemp[ this->modChannelIndices[i] ].hhDim;
    this->chanModDims[i].phDim = chanModDimsTemp[ this->modChannelIndices[i] ].phDim;
    this->chanModDims[i].hpDim = chanModDimsTemp[ this->modChannelIndices[i] ].hpDim;
  }

  for(std::size_t i = 0; i < this->nSpstates; i++){
    this->threeBodyChanDims[i].phhDim = threeBodyChanDimsTemp[ this->spChannelIndices[i] ].phhDim;
    this->threeBodyChanDims[i].hppDim = threeBodyChanDimsTemp[ this->spChannelIndices[i] ].hppDim;
  }

  // printf("CHECK\n");
  // for(std::size_t i = 0; i < this->nChannels; i++){
  //   printf("i: %zu, reg of inv: %zu\n",i, channelIndices[ inverseChannelIndices[i] ]);
  // }

  delete [] chanDimsTemp;
  delete [] chanModDimsTemp;
  delete [] threeBodyChanDimsTemp;
}
