DATA_SECTION
//Pribilof Island red king crab model

  !!CLASS ofstream post("eval.csv");
  int call_no;
  !! call_no = 0;
   int endyrTrawl; 
   int Nsex;
  !! Nsex = 2;                                                     		 // Number of sex

  init_int styr                                                    			 // start year of the model
  init_int endyr															//end year of the model
 !! endyrTrawl = endyr-1;
 
  // Survey data
  init_int lenbinN                                                    		// number of length bins for males in the model
  init_int nobs_survey                                                		// number of years of survey numbers data
  init_ivector yrs_survey(1,nobs_survey)                            // years when have survey data
  init_int nobs_surv_len                                              		// number of years of survey length data
  init_ivector yrs_survey_len(1,nobs_survey)                     // years when have fsurvey length data
  init_ivector nsamples_surv_M(1,nobs_survey)   	            // nsamples weight for survey length comps M
  init_ivector nsamples_surv_F(1,nobs_survey)   	            // nsamples weight for survey length comps F
  
  init_vector SurvMnum(1,nobs_survey)							//survey index of abundance (numbers) males
  init_vector SurvMcv(1,nobs_survey)								//CV for male survey index
  init_vector SurvFnum(1,nobs_survey)							    //survey index of abundance (numbers) females
  init_vector SurvFcv(1,nobs_survey)								//CV for female survey index
  !!cout<<SurvFcv<<"SurvF"<<endl;
  init_matrix SurvMlen(1,nobs_survey,1,lenbinN)					//male survey length frequencies
  init_matrix SurvFlen(1,nobs_survey,1,lenbinN)					//female survey length frequencies    
  init_vector length_bins(1,lenbinN)									//length bins
  init_matrix PropMature(1,2,1,lenbinN)								//proportion mature males
 !!cout<<"PropMat"<<PropMature<<endl;
    // Catch data
  !! ad_comm::change_datafile_name("RKC Catches.DAT");
  init_int nobs_catch														//number of years for which catch data are available
  init_ivector yrs_catch(1,nobs_catch)								//years for which catch data are available
  init_vector CatchNum(1,nobs_catch)								//catch in numbers
  init_vector CatchBio(1,nobs_catch)								//catch in biomass (t)

  init_int nobs_disc															//number of years for which discard data are available
  init_ivector yrs_disc(1,nobs_disc)									//years for which discard data are availabale
  init_vector discLegalM(1,nobs_disc)								//discard of legal males in biomass (t)
  init_vector discSubM(1,nobs_disc)								//discard of sublegal males in biomas (t)
  init_vector discF(1,nobs_disc)										//discard of females in biomass (t)
 
  init_int nobs_bycatch													//number of years for which bycatch data are available 
  init_ivector yrs_bycatch(1,nobs_bycatch)						//years for which bycatch data are available
  init_vector FixedGearBycatch(1,nobs_bycatch)				//fixed gear bycatch biomass in (t)
  init_vector TrawlBycatch(1,nobs_bycatch)						//trawl bycatch in (t)
  init_vector TrawlBycatch20off(1,nobs_bycatch)						//trawl bycatch in (t)  
   !!cout<<"TrawlBycatch"<<TrawlBycatch20off<<endl;
  
    // Control file
  !! ad_comm::change_datafile_name("PIRKCv2.ctl");
  init_number smallNum
  init_number natM
  init_number molt_slope
  init_number molt_L50  
  init_number Catcv
  init_number byCatcv
  init_number InitSmoothWeight
  init_number ShrinkCV
  init_number weightAm
  init_number weightBm
  init_number weightAf
  init_number weightBf
  init_number trawl_sel50
  init_number trawlSlope
  init_number trawlPenSmooth
  init_number FitFemales
  init_number TimeVaryQ
  init_number TimeVaryQsmooth
  init_number AdjustEffN
  init_number BBRKC_selUse
  init_number BBRKC_sel50
  init_number BBRKC_slope
  init_number BBRKC_discMort
  init_number BBRKC_discOmega
  init_number BBRKC_discKappa
  init_number BBRKC_discGamma

  int ipass;
  int Nproj;
  int IsB0;                                                         // Set to 0 for B0
  !! Nproj = 125;

// =======================================================================
PARAMETER_SECTION

    init_bounded_vector srv_q(styr,endyr,0.1,1,3)
	init_bounded_number fish_sel50(137.9,138.1,-2)
	init_bounded_number fish_sel95(137.95,138.15,-2)
	init_bounded_number srv_sel50(50,250,1)
	init_bounded_number srv_sel95(60,250,1)

 //estimated Fs
   init_number log_avg_fmort_dir(1)
   init_bounded_dev_vector fmort_dir_dev(1993,1998,-100,100,1)
   init_number log_avg_fmort_trawl(1)
   init_bounded_dev_vector fmort_trawl_dev(1991,endyrTrawl,-100,100,1)
   // init_bounded_vector stNatLenM(1,lenbinN,-100,100,1)
   // init_bounded_vector stNatLenF(1,lenbinN,-100,100,1)
   
   init_number mean_log_rec(1)
   init_bounded_dev_vector rec_dev(styr,endyr,-200,200,1)
   
   init_bounded_number af(0,100,2)
   init_bounded_number am(0,100,2)
   init_bounded_number bf(.7,1.5,2)
   init_bounded_number bm(.7,1.5,2)
   init_bounded_vector growth_beta(1,2,0.1,3,3)
   init_bounded_number alpha1_rec(.5,10,3)
   init_bounded_number beta_rec(.1,9,3)
   // init_bounded_number dummy(1,2)
  
  //end of estimated parameters///
  objective_function_value f
  
  //Biological processes (first index is males, second is females)
   matrix molt_p(1,2,1,lenbinN)											//molting probabilities
   matrix mean_length(1,2,1,lenbinN)                                // Predicted post-moult sizes
   vector rec_len(1,lenbinN)                                          	 // Recruitment length
  3darray len_len(1,2,1,lenbinN,1,lenbinN)                         // length to length growth array
 
  vector sel_fish(1,lenbinN)											    // retained selectivity
  vector sel_fish_discard(1,lenbinN)								    // fishery discard selectivity
  matrix sel_srv(styr,endyr+Nproj,1,lenbinN)            	    	// Survey selectivity (do I need different eras?)
  vector sel_trawl(1,lenbinN)											// selectivity for trawl bycatch
   
  //fishing mortality and survival by area
  vector fmort_dir(styr,endyr+Nproj)								//Directed fully selected F on males
  vector fmort_trawl(styr,endyr+Nproj)						    //Trawl fully selected F on males
  matrix F_dir(styr,endyr+Nproj,1,lenbinN)
  matrix F_disc(styr,endyr+Nproj,1,lenbinN)  
  matrix F_trawl(styr,endyr+Nproj,1,lenbinN)
   
  //Numbers at length by sex (ladies first)
  3darray Natlen(1,2,styr,endyr+Nproj,1,lenbinN)
  sdreport_matrix Spbio(1,2,styr,endyr+Nproj)

 //predicted total numbers in the survey
   vector predSurvM(styr,endyr+Nproj)
   vector predSurvF(styr,endyr+Nproj)
   matrix predSurvNatlenM(styr,endyr+Nproj,1,lenbinN)
   matrix predSurvNatlenF(styr,endyr+Nproj,1,lenbinN)
   
  //predicted catches, discards and bycatches at length
  3darray predCatchAtLen(1,2,styr,endyr+Nproj,1,lenbinN)
  3darray PredTrawlAtLen(1,2,styr,endyr+Nproj,1,lenbinN)
  
  //pred catch N (for aggregating CatLen)
  sdreport_vector predCatchN(1,nobs_catch);  
  sdreport_matrix predTrawlCatchN(1,2,styr,endyr+Nproj);  
  
  matrix predSurvLenFreqM(styr,endyr+Nproj,1,lenbinN)
  matrix predSurvLenFreqF(styr,endyr+Nproj,1,lenbinN)

  //likelihood components
  number SurvM_like
  number SurvF_like
  number SurvMLen_like
  number SurvFLen_like
  number Catch_like
  number Trawl_like
  number initsmo_penalF
  number initsmo_penalM
  number Rec_pen
  number F_pen
  number trawlF_pen
  number VaryQ_pen
 
  //weight at length
  matrix WeightAtLen(1,2,1,lenbinN)
  
  //reference points
  number F35;                                                       // F35
  number SBPRF35;                                                   // SBPR(F=F35);
  number Bmsy;                                                      // Bmsy
  number FOFL;														//FOFL
  number OFL;
  number FutRec;                                                     // Future recruitment
  number FutMort;                                                    // Future target
  number FOFL_tier4;														//FOFL
  number OFL_tier4;
  number Bmsy_tier4; 

//  ========================================================================
PRELIMINARY_CALCS_SECTION
 int i;
 
 for(i=1;i<=lenbinN;i++)
 {
  WeightAtLen(2,i) = weightAm * pow(length_bins(i),weightBm);
  WeightAtLen(1,i) = weightAf * pow(length_bins(i),weightBf);
  }

  //  ========================================================================
INITIALIZATION_SECTION
 af 3
 am 3.5
 bf 1.05
 bm 1.15
 growth_beta 0.75
 alpha1_rec 4
 beta_rec .4
// ==========================================================================
PROCEDURE_SECTION
   getmoltingp();
   // cout<<"Molt"<<endl;
  getgrowth();
   // cout<<"Growth"<<endl;
  getselectivity();
   // cout<<"Select"<<endl;
  getmortality();
  // cout<<"Mort"<<endl;
  get_num_at_len();
  // cout<<"Numbers"<<endl;
  evaluate_the_objective_function();

    if (mceval_phase())
   {
    Find_F35();
    Find_OFL();
    Find_OFL_tier4();
	
    post << Bmsy << " " << F35 << " " << FOFL << " " << OFL << " ";   
    post << srv_sel50 << " " << srv_sel95 << " " << log_avg_fmort_dir << " " << fmort_dir_dev << " " << log_avg_fmort_trawl << " ";
    post <<  fmort_trawl_dev << " " << mean_log_rec << " " << rec_dev << " ";
	post << af << " " <<am << " "<< bf << " " << bm <<" ";
	post << growth_beta << " " <<alpha1_rec << " "<< beta_rec << " " ;
	post << Bmsy_tier4 << " " <<FOFL_tier4 << " "<< OFL_tier4 << " "  <<f <<" " <<endl;
   }
  
// ==========================================================================
FUNCTION getmoltingp
 int j;
 //females molt annually
  molt_p(1)= 1;
  //males molt annually until maturity at which point molting become less frequent
  for(j=1;j<=lenbinN;j++)
   molt_p(2,j) = 1 - 1./(1.+mfexp(-1.*molt_slope*(length_bins(j)-molt_L50)));

// ==========================================================================
FUNCTION getgrowth
 int ilen,il2,sex;
  dvariable devia, recsum, alpha_rec,alpha,lensum;

  //mean_length is the expected size after molting from a given size bin
  for(ilen=1;ilen<=lenbinN;ilen++)
   { 
      // mean_length(1,ilen)= length_bins(ilen) + grow_Q_f/(1+mfexp(-log(19)*(length_bins(ilen)-grow_l50_f)/(grow_l95_f-grow_l50_f)));
      // mean_length(2,ilen)= am + bm * length_bins(ilen);
      mean_length(1,ilen)= af + bf * length_bins(ilen);
      mean_length(2,ilen)= am + bm * length_bins(ilen);
   }
   
   // cout<<"meanLen"<<mean_length<<endl;
// using Gamma function for transition matrix
// devia is the bounds of growth bins to evaluate
// the gamma function (x) in prop = integral(i1 to i2) g(x|alpha,beta) dx
// alpha and growth_beta are parameters 
// alpha is the mean growth increment per molt for some premolt length class
// alpha = mean growth increment per molt divided by beta
// beta is the shape parameter - larger beta - more variance 
 
  for (sex=1;sex<=2;sex++)
   for(ilen=1;ilen<=lenbinN;ilen++)
    {

     // subract the 2.5 from the midpoint of the length bin to get the lower bound
     alpha = (mean_length(sex,ilen)-(length_bins(ilen)-2.5))/growth_beta(sex);
    
     lensum = 0;
     for(il2=1;il2<=lenbinN;il2++)
      if(il2>=ilen)
       {
        devia = length_bins(il2)+2.5-length_bins(ilen);
        len_len(sex,ilen,il2) = mfexp((alpha-1.0)*log(devia)-devia/growth_beta(sex));
        lensum += len_len(sex,ilen,il2);
       }  

    //standardize so each row sums to 1.0
    for(il2=1;il2<=lenbinN;il2++)
      len_len(sex,ilen,il2)=len_len(sex,ilen,il2)/lensum;
   }
   
 // Fraction recruiting
  recsum=0.0;
  alpha_rec=alpha1_rec/beta_rec;
  for(ilen=1;ilen<=lenbinN;ilen++)
   {
    devia = length_bins(ilen)+2.5-length_bins(1);
    rec_len(ilen) = pow(devia,alpha_rec-1.)*mfexp(-devia/beta_rec);
    recsum += rec_len(ilen);
   }

  //standardize so each row sums to 1.0
  for(ilen=1;ilen<=lenbinN;ilen++)
   rec_len(ilen) = rec_len(ilen)/recsum;


   
// ==========================================================================
FUNCTION getselectivity
 int j,i;
 dvariable tempSel;

  // trawl selectivity curves
    for (j=1;j<=lenbinN;j++)
     sel_trawl(j)     = 1./(1.+mfexp(-1.*trawlSlope* (length_bins(j)-trawl_sel50))); 

 
  // fishery selectivity curves (knife-edge at legal size, no discard mort)
  if(BBRKC_selUse==0)
  {
    for (j=1;j<=lenbinN;j++)
    sel_fish(j)       =  1./(1.+mfexp(-1.*log(19) *  (length_bins(j)-fish_sel50)/(fish_sel95-fish_sel50))); 
   }

  // fishery selectivity curves
  if(BBRKC_selUse>0)
  {
    for (j=1;j<=lenbinN;j++)
	{
    sel_fish(j)       		=  1./(1.+mfexp(-1*BBRKC_slope*(length_bins(j)-BBRKC_sel50))); 
	
	// this is silly
	if(lenbinN<=135)
	{
	 tempSel = BBRKC_discOmega + BBRKC_discKappa*length_bins(j);
	  sel_fish_discard(j)  =0;
     if(tempSel>0)
      sel_fish_discard(j)  =  tempSel;
    }
	if(lenbinN>135)
      sel_fish_discard(j)  =  sel_fish_discard(j-1) + 5 * BBRKC_discGamma;
	}
  }
    
 //survey selectivity curves 
  for(i=styr;i<=endyr;i++)	 
   for(j=1;j<=lenbinN;j++)
   {
	if(TimeVaryQ>0)
     sel_srv(i,j)		   = srv_q (i)/ (1.+mfexp(-1.*log(19.)*(length_bins(j)-srv_sel50)/(srv_sel95-srv_sel50)));
	if(TimeVaryQ==0)
     sel_srv(i,j)		   = 1 / (1.+mfexp(-1.*log(19.)*(length_bins(j)-srv_sel50)/(srv_sel95-srv_sel50)));
    }
	
  for(i=endyr+1;i<=endyr+Nproj;i++)	 
   for(j=1;j<=lenbinN;j++)
   {
	if(TimeVaryQ>0)
     sel_srv(i,j)		   = srv_q (endyr)/ (1.+mfexp(-1.*log(19.)*(length_bins(j)-srv_sel50)/(srv_sel95-srv_sel50)));
	if(TimeVaryQ==0)
     sel_srv(i,j)		   = 1 / (1.+mfexp(-1.*log(19.)*(length_bins(j)-srv_sel50)/(srv_sel95-srv_sel50)));
    }
  // ==========================================================================
FUNCTION getmortality
 int i,j;
 //Fourier forgive my hardcoding sins

 //directed fishery retained and discard mortality
 for (i=1993;i<=1998;i++)
 {
   fmort_dir(i) 		= mfexp(log_avg_fmort_dir    + fmort_dir_dev(i));
   F_dir(i)			= sel_fish*fmort_dir(i);
   F_disc(i) 		= sel_fish_discard*BBRKC_discMort;
  }

 for(i=styr;i<=1992;i++)
 {
  F_dir(i) = 0; 
  F_disc(i) = 0; 
 }

 for(i=1999;i<=endyr;i++)
 {
  F_dir(i) = 0; 
  F_disc(i) = 0; 
 }
 
 for (i=styr;i<=endyr;i++)
    F_trawl(i) = 0; 	
   
 for (i=1991;i<=endyrTrawl;i++)
  {
    fmort_trawl(i)  = mfexp(log_avg_fmort_trawl + fmort_trawl_dev(i));
    F_trawl(i) = sel_trawl*fmort_trawl(i); 	
   } 
   

   
// ==========================================================================
FUNCTION get_num_at_len
//to estimate an initial recruitment in the start year?
 int i;
 Natlen.initialize();
 Spbio.initialize();
 predSurvM.initialize();
 predSurvF.initialize();
 
 // Natlen(1,styr) = mfexp(stNatLenF);
 // Natlen(2,styr)  = mfexp(stNatLenM);

 // Natlen(1,styr) = mfexp(mean_log_rec)*rec_len;
 // Natlen(2,styr)  = mfexp(mean_log_rec)*rec_len;  

 Natlen(1,styr) = mfexp(mean_log_rec+rec_dev(styr))*rec_len;
 Natlen(2,styr)  = mfexp(mean_log_rec+rec_dev(styr))*rec_len; 
 
  for(ipass=styr;ipass<=endyr;ipass++) 
    get_num_at_len_yr();

	// ==========================================================================
FUNCTION get_num_at_len_yr
 int i,j,sex;
 dvar_vector Molting(1,lenbinN);
 dvar_vector NotMolting(1,lenbinN);
 dvar_vector tempMolt(1,lenbinN);
 dvar_vector tempN(1,lenbinN);
 dvar_vector tempNmate(1,lenbinN);
 dvar_vector tempBio(1,lenbinN);
 dvariable surv2fish,surv2surv,surv2mate;
 
 i = ipass;
  // cout<<"NTLEN"<<Natlen(1,i)<<endl;
  // cout<<"NTLEN"<<Natlen(2,i)<<endl; 

  //survey occurs in July
  predSurvNatlenF(i) = elem_prod(sel_srv(i),Natlen(1,i));
  predSurvNatlenM(i)  = elem_prod(sel_srv(i),Natlen(2,i));
  
  for(j=1;j<=lenbinN;j++)
   {
     predSurvM(i) += predSurvNatlenM(i,j);
     predSurvF(i)  += predSurvNatlenF(i,j);
   }	 
  // cout<<"PredSurvM"<<predSurvM(i)<<endl;
   //survey len freqs
  predSurvLenFreqM(i) = predSurvNatlenM(i)/predSurvM(i);
  predSurvLenFreqF(i) = predSurvNatlenF(i)/predSurvF(i);
 if(i<=endyr)  
  surv2fish = mfexp(-3*natM/12) ;		//fishery historically happened 3 months after survey (June 15 to September 15)
 if(i>endyr)  
  surv2fish = mfexp(-4*natM/12) ;		//fishery now happens 4 month after survey (June 15 to October 15)
  surv2surv = mfexp(-9*natM/12);	//survey happens 9 months after fishery 
  surv2mate = mfexp(-5*natM/12);	//mating happens 5 months after fishing
 // cout<<"surv2surv"<<surv2surv<<endl;
 for(sex=1;sex<=2;sex++) 
  {
   //ladies first
    tempN = surv2fish * Natlen(sex,i);
   // cout<<"tempN"<<tempN<<endl;
    predCatchAtLen(sex,i) = elem_prod(1-mfexp(-F_dir(i)),tempN);								//predicted catch at length

    tempN = elem_prod(mfexp(-F_dir(i)),tempN);														//fishery occurs
	   // cout<<"tempN"<<tempN<<endl;
    tempN  = elem_prod(mfexp(-F_disc(i)),tempN);														//fishery discard mortality occurs
	   // cout<<"tempN3"<<tempN<<endl;
    tempN  = elem_prod(mfexp(-F_trawl(i)),tempN);													//bycatch occurs
	   // cout<<"tempN3"<<tempN<<endl;
    PredTrawlAtLen(sex,i) = elem_prod(1-mfexp(-F_trawl(i)),tempN);							//predicted bycatch at length

	tempBio=0;																										//predicted biomass of bycatch (t)
	tempBio= elem_prod(PredTrawlAtLen(sex,i),WeightAtLen(sex));
    predTrawlCatchN(sex,i) = 0;
   for(j=1;j<=lenbinN;j++)
	predTrawlCatchN(sex,i) +=tempBio(j)/1000000;

   // cout<<"Predtrawl"<<predTrawlCatchN(sex,i)<<endl;	
   Molting = elem_prod(molt_p(sex),tempN);																//molting
   NotMolting = elem_prod(1-molt_p(sex),tempN);
   tempMolt = Molting * len_len(sex);

   tempN = NotMolting + tempMolt;
   	   // cout<<"tempN4"<<tempN<<endl;
   if (i <=endyr)																									   //recruitment
    tempN += mfexp(mean_log_rec + rec_dev(i))*rec_len;
   else 
    tempN += FutRec*rec_len;
	   // cout<<"tempN5"<<tempN<<endl;
    Spbio(sex,i) = 0;		//calculate spawning biomass
	tempNmate = tempN*surv2mate;
    tempBio = elem_prod(tempNmate,PropMature(sex));
	tempBio = elem_prod(tempBio,WeightAtLen(sex));
	
   for(j=1;j<=lenbinN;j++)
    Spbio(sex,i) += tempBio(j);

	// cout<<"spbio"<<Spbio(sex,i)<<endl;
   Natlen(sex,i+1) = surv2surv * tempN;																	//natural mortality to survey
   }

 //fishery occurs in September (a month after the survey)
 //take out discard with fishery (this will influence the magnitude of the F associated
 // discards and bycatch, but it will not influence calculation of management quantities) 
 //growth occurs in the winter
 //recruitment occurs before survey
 //mating occurs in February (7 months after the survey)
// ==========================================================================                                                             
FUNCTION evaluate_the_objective_function
 int i,j;
 
 SurvM_like.initialize();
 SurvF_like.initialize();
 SurvMLen_like.initialize();
 SurvFLen_like.initialize();
 Catch_like.initialize();
 Trawl_like.initialize();
 initsmo_penalM.initialize();
 initsmo_penalF.initialize();
 trawlF_pen.initialize();
 F_pen.initialize();

 //Survey abundance
	for(i=2;i<=nobs_survey;i++)
	 {
	 SurvM_like +=  square(log(SurvMnum(i) + smallNum) - log(predSurvM(yrs_survey(i)) + smallNum))  / ((2) * (log(ShrinkCV*SurvMcv(i)*SurvMcv(i) + 1)  ));
	 SurvF_like  +=  square(log(SurvFnum(i) + smallNum) - log(predSurvF(yrs_survey(i)) + smallNum))   / ((2) * (log(ShrinkCV*SurvFcv(i)*SurvFcv(i) + 1)  ));
	 }

   //survey length frequencies
	for(i=1;i<=nobs_survey;i++)
	 for(j=1;j<=lenbinN;j++)
	 {
	 if(SurvMlen(i,j)>0.001)
	  SurvMLen_like -= AdjustEffN*nsamples_surv_M(i) *SurvMlen(i,j)   *log(predSurvLenFreqM(yrs_survey(i),j) );
	 if(SurvFlen(i,j)>0.001)
	  SurvFLen_like  -= AdjustEffN*nsamples_surv_F(i) *SurvFlen(i,j)    *log(predSurvLenFreqF(yrs_survey(i),j) );
	  }
 //catch numbers
  predCatchN.initialize();
  for(i=1;i<=nobs_catch;i++)
   for(j=1;j<=lenbinN;j++)
    predCatchN(i) += predCatchAtLen(2,yrs_catch(i),j);

	for(i=1;i<=nobs_catch;i++)
  	Catch_like +=  square(log(CatchNum(i) + smallNum) - log(predCatchN(i) + smallNum))  / (sqrt(2) * sqrt(log(Catcv*Catcv+ 1)  ));
 
    for(i=1;i<=nobs_bycatch;i++)
    	Trawl_like +=  square(log(TrawlBycatch20off(i) + smallNum) - log(predTrawlCatchN(1,yrs_bycatch(i)) + predTrawlCatchN(2,yrs_bycatch(i)) + smallNum))  / (sqrt(2) * sqrt(log(byCatcv*byCatcv+ 1)  ));

  Rec_pen = InitSmoothWeight*(norm2(first_difference(rec_dev)));
  F_pen = InitSmoothWeight*(norm2(first_difference(fmort_dir_dev)));
  trawlF_pen = trawlPenSmooth*(norm2(first_difference(fmort_trawl_dev)));
  VaryQ_pen = TimeVaryQsmooth*(norm2(first_difference(srv_q(styr,endyr))));

  
 f = 0;
 f += SurvM_like;
 f += SurvMLen_like;
 if(FitFemales>0)
 {
  f += SurvF_like;
  f += SurvFLen_like;
 }
 f += Catch_like; 
 f += Trawl_like;
 f += Rec_pen;
 f += F_pen;
 f += trawlF_pen;
 if(TimeVaryQ>0)
 f += VaryQ_pen;
 cout << current_phase() << " " << call_no << " " << f << endl;
 
 // ==========================================================================

FUNCTION get_fut_mortality
 int i,j,nn;
 dvariable AvgTrawlF;
 
 nn = 0;
 AvgTrawlF =0;
 for(j=(endyr-6);j<endyr;j++)
  {
    AvgTrawlF += mfexp(log_avg_fmort_trawl + fmort_trawl_dev(j));
    nn += 1;
  }
  AvgTrawlF = AvgTrawlF/nn;
  
  for (i=ipass;i<=endyr+100;i++)
  {
   if (IsB0 == 0)
    {
   	fmort_dir(i) = 0; 
	fmort_trawl(i)  = 0;
	}
   else 
    {
	fmort_dir(i) = FutMort; 
    fmort_trawl(i)  = AvgTrawlF;
	}
    F_trawl(i) = sel_trawl*fmort_trawl(i); 	
    F_dir(i) = sel_fish*fmort_dir(i);
   }
 
// ==========================================================================

FUNCTION Find_F35
  dvariable Ratio,Target,B0;
  dvariable Ntot,totB0,Btest;
  int icnt,i;
  dvar_vector tempN(1,lenbinN);

  // Find B0
  IsB0 = 0;
  FutRec = 100000;
  FutMort = 0;
  ipass = endyr+1;
  get_fut_mortality();
  for (ipass=endyr+1;ipass<=endyr+100;ipass++) get_num_at_len_yr(); 
  
 // put biomass in here soon
  B0 = Spbio(2,endyr +99);
  cout<<"B0"<<B0<<endl;
  Target = 0.35;
  
// Find F35%  
  IsB0 = 1;
  FutMort = 0.3;
  for (icnt=1;icnt<=20;icnt++)
   {
    ipass = endyr+1;
    get_fut_mortality();
    for (ipass=endyr+1;ipass<=endyr+100;ipass++) get_num_at_len_yr(); 
     Btest = Spbio(2,endyr +99);
   
   Ratio = Btest/B0;
    cout << FutMort << " " << Ratio << endl;
    FutMort = FutMort * Ratio / Target;
   }
  F35 = FutMort; 
  SBPRF35 = Btest/FutRec;
   cout << FutMort << " " << Ratio << endl; 

// ==========================================================================

FUNCTION Find_OFL
  dvariable Fmsy,Rbar,nn,alpha,beta;
  int BMSY_Yr1, BMSY_Yr2,ii,Iyr,kk,jj;
  
  BMSY_Yr1 = 1983;BMSY_Yr2 = endyr-1;
  alpha = 0.05;
  beta = 0.25;
  
  // Define Fmsy
  Fmsy = F35;
  
 // Find Rbar (Dynamic or not)
  Rbar = 0; nn= 0;
  for (Iyr=BMSY_Yr1;Iyr<=BMSY_Yr2;Iyr++)
   {
    Rbar += mfexp(mean_log_rec+rec_dev(Iyr));
    nn += 1;
   }
  Rbar = Rbar / nn;

 // Specify the BMSY proxy
  Bmsy = SBPRF35 * Rbar;
  cout << "Rbar" << Rbar << endl;
  cout << "Bmsy" << Bmsy << endl;

  cout<<"FMSY"<<Fmsy<<endl;
  // year for projection
  ipass = endyr+1;

  // Define future recruitment 
  if (ipass > endyr) FutRec = Rbar;
 
  // Find FOFL
  FutMort = Fmsy;
  get_fut_mortality();
  get_num_at_len_yr();
      cout<<"Mspbio"<<Spbio(2,ipass)<<endl;
  if (Spbio(2,ipass) < Bmsy)
   {
    FutMort = 0;
    get_fut_mortality();
    get_num_at_len_yr();
    if (Spbio(2,ipass) > beta*Bmsy)
     {
      FutMort = Fmsy/2;
      get_fut_mortality();
      get_num_at_len_yr();
      for (ii=1;ii<=10;ii++)
       {
        FutMort = Fmsy*(Spbio(2,ipass)/Bmsy-alpha)/(1-alpha);
        get_fut_mortality();
        get_num_at_len_yr();
       }
     }
   }

   FOFL = FutMort;
   cout<<"FOFL"<<FOFL<<endl;
   get_fut_mortality();
   get_num_at_len_yr();
   OFL = 0;
   for(ii = 1;ii<=lenbinN;ii++)
    OFL += predCatchAtLen(2,ipass,ii) *WeightAtLen(2,ii);
   cout<<"OFL"<<OFL<<endl;
// ==========================================================

FUNCTION Find_OFL_tier4
  dvariable Fmsy,Rbar,nn,alpha,beta;
  int BMSY_Yr1, BMSY_Yr2,ii,Iyr,kk,jj;
  
 BMSY_Yr1 = 1991;BMSY_Yr2 = endyr-1;
  alpha = 0.05;
  beta = 0.25;
  
  // Define Fmsy tier 4
  Fmsy = natM;
  
   // Find Rbar (Dynamic or not)
  Rbar = 0; nn= 0;
  for (Iyr=BMSY_Yr1;Iyr<=BMSY_Yr2;Iyr++)
   {
    Rbar += mfexp(mean_log_rec+rec_dev(Iyr));
    nn += 1;
   }
  Rbar = Rbar / nn;
  
 // Specify the BMSY proxy (tier 4)
  Bmsy_tier4 = 0; nn= 0;
  for (Iyr=BMSY_Yr1;Iyr<=BMSY_Yr2;Iyr++)
   {
    Bmsy_tier4 += Spbio(2,Iyr);
    nn += 1;
   }
  Bmsy_tier4 = Bmsy_tier4 / nn;

  // year for projection
  ipass = endyr+1;

  // Define future recruitment 
  if (ipass > endyr) FutRec = Rbar;
 
  // Find FOFL
  FutMort = Fmsy;
  get_fut_mortality();
  get_num_at_len_yr();
  if (Spbio(2,ipass) < Bmsy)
   {
    FutMort = 0;
    get_fut_mortality();
    get_num_at_len_yr();
    if (Spbio(2,ipass) > beta*Bmsy)
     {
      FutMort = Fmsy/2;
      get_fut_mortality();
      get_num_at_len_yr();
      for (ii=1;ii<=10;ii++)
       {
        FutMort = Fmsy*(Spbio(2,ipass)/Bmsy-alpha)/(1-alpha);
        get_fut_mortality();
        get_num_at_len_yr();
       }
     }
   }

   FOFL_tier4 = FutMort;
   cout<<"FOFL"<<FOFL<<endl;
   get_fut_mortality();
   get_num_at_len_yr();
   OFL_tier4 = 0;
   for(ii = 1;ii<=lenbinN;ii++)
    OFL_tier4 += predCatchAtLen(2,ipass,ii) *WeightAtLen(2,ii);
   cout<<"OFL"<<OFL_tier4<<endl;
// ==========================================================


REPORT_SECTION
 int i,j;
 
  Find_F35();
  Find_OFL();
  Find_OFL_tier4();
  
  report<<"SBPRF35"<<endl;
  report<<SBPRF35<<endl; 
  report<<"B35"<<endl;
  report<<Bmsy<<endl;
  report<<"F35"<<endl;
  report<<F35<<endl;
  report<<"FOFL"<<endl;
  report<<FOFL<<endl;
  report<<"OFL"<<endl;
  report<<OFL<<endl;
  report<<"Bmsy_tier4"<<endl;
  report<<Bmsy_tier4<<endl;
  report<<"FOFL_tier4"<<endl;
  report<<FOFL_tier4<<endl;
  report<<"OFL_tier4"<<endl;
  report<<OFL_tier4<<endl;
  
 report << "#pred male survey N"<<endl;
 report << predSurvM << endl;
 report << "#pred female survey N"<<endl;
 report << predSurvF << endl;
 report << "#obs male survey N"<<endl;
 report << SurvMnum << endl;
 report << "#obs female survey N"<<endl;
 report << SurvFnum << endl;
 
 report << "#pred catch N"<<endl;
 report << predCatchN << endl;
 report << "#obs catch N"<<endl;
 report << CatchNum << endl;
 
 report << "#obs surv M at len" <<endl;
    for(j=1;j<=nobs_survey;j++)
     report <<SurvMlen(j) << endl; 
 report << "#obs surv F at len" <<endl;
    for(j=1;j<=nobs_survey;j++)
     report <<SurvFlen(j) << endl;  
	 
 report << "#pred surv M at len" <<endl;
    for(j=styr;j<=endyr;j++)
     report <<predSurvLenFreqM(j) << endl; 
 report << "#pred surv F at len" <<endl;
    for(j=styr;j<=endyr;j++)
     report <<predSurvLenFreqF(j) << endl;  	 
 
  report<<"#male spawning biomass"<< endl;
  report<<Spbio(2)<<endl;
  report<<"#female spawning biomass"<< endl;
  report<<Spbio(1)<<endl;
  report<<"#obs trawl bycatch"<< endl; 
   report<<TrawlBycatch20off<<endl;
  report<<"#pred trawl bycatch fem"<< endl;  
  report<<predTrawlCatchN(1)<<endl;
  report<<"#pred trawl bycatch male"<< endl;  
  report<<predTrawlCatchN(2)<<endl;
  
  report<<"#size transition matrix Males"<<endl;
   for(j=1;j<=lenbinN;j++)
    report<<len_len(2,j)<<endl;
	
  report<<"#fraction recruiting"<<endl;
  report<<rec_len<<endl;
 
   report<<"#size transition matrix females"<<endl;
   for(j=1;j<=lenbinN;j++)
    report<<len_len(1,j)<<endl;
 
//===============================================================================

GLOBALS_SECTION
 #include <math.h>
 #include <admodel.h>
  #include <time.h>
 
 ofstream CheckFile;
  time_t start,finish;
  long hour,minute,second;
  double elapsed_time;

// EVALUATE THE GAMMA FUNCTION
 dvariable gammln(dvariable& xx)
 {
   RETURN_ARRAYS_INCREMENT();	
        dvariable x,y,tmp,ser;
        static double cof[6]={76.18009172947146,-86.50532032941677,
                24.01409824083091,-1.231739572450155,
               0.1208650973866179e-2,-0.5395239384953e-5}; 
        int j;
        y=x=xx;
        tmp=x+5.5;
        tmp -= (x+0.5)*log(tmp);
        ser=1.000000000190015;
        for (j=0;j<=5;j++) 
          {
            ser += cof[j]/(y+1.);
          }
         dvariable value_=-tmp+log(2.5066282746310005*ser/x);
    RETURN_ARRAYS_DECREMENT();         
         return(value_); 
 }

// ===============================================================================
RUNTIME_SECTION
//one number for each phase, if more phases then uses the last number
  maximum_function_evaluations 5000,5000,5000,5000,5000,5000,5000
//  convergence_criteria 1,1,1,1,.01,.001,1e-3,1e-3
 convergence_criteria .1,.001,.001

TOP_OF_MAIN_SECTION
  arrmblsize = 4000000;
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(3000000); // this may be incorrect in
  // the AUTODIF manual.
  gradient_structure::set_CMPDIF_BUFFER_SIZE(100000000);
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(900);
  
  time(&start);
  CheckFile.open("Check.Out");

FINAL_SECTION
 time(&finish); 
 elapsed_time = difftime(finish,start);
 hour = long(elapsed_time)/3600;
 minute = long(elapsed_time)%3600/60;
 second = (long(elapsed_time)%3600)%60;
 
 cout<<"Survey M"<<SurvM_like<<endl;
 cout<<"Survey F"<<SurvF_like<<endl;
 cout<<"Survey M len"<<SurvMLen_like<<endl;
 cout<<"Survey F len"<<SurvFLen_like<<endl;
 cout<<"Catch"<<Catch_like<<endl;
 cout<<"Trawl"<<Trawl_like<<endl;
 cout<<"VaryQ_pen"<<VaryQ_pen<<endl; 
 cout<<"trawlF_pen"<<trawlF_pen<<endl; 
 cout<<"F_pen"<<F_pen<<endl; 
 cout<<"Rec_pen"<<Rec_pen<<endl; 

 
 cout << endl << endl << "Starting time: " << ctime(&start);
 cout << "Finishing time: " << ctime(&finish);
 cout << "This run took: " << hour << " hours, " << minute << " minutes, " << second << " seconds." << endl << endl;
