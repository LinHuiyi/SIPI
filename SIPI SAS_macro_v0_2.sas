/*****************
•	SAS macro name	: SIPI 
•	Version		: v0.2 
•	Contents		: SNP Interaction Pattern Identifier (SIPI) for the binary outcome
•	Date: 9/21/2016

SIPI is used to evaluate pairwise SNP-SNP interactions associated with a binary outcome. For each SNP pair, the SIPI evaluates 45 logistic interaction models which take inheritance modes [additive, dominant, and recessive mode with the original (based on the minor allele) and reverse coding], and risk category grouping (model structure: hierarchical and non-hierarchical models) into consideration. The best interaction pattern is the one with the lowest value of the Bayesian information criterion (BIC). The significant test of the interaction is based on the Wald test (or the likelihood ratio test) of the interaction term in a model. There are two sub-functions of this SIPI macro: (1) one-pair only and (2) pairwise interaction analyses. The details of the 45 models/patterns are listed in the SIPI paper.
There are two sub-functions of this SIPI macro:
(1) For testing interaction for only one SNP pair, the outputs include (a) SNP description statistics table and (b) a list of all 45 model results. 
(2) For pairwise SNP interactions, the outputs include (a) SNP description statistics table and (b) a list of the best interaction pattern for each SNP pair.

•	Macro developers: 
Hui-Yi Lin, PhD
email: hlin1@lsuhsc.edu

Yung-Hsin Liu, MS
email: James.Liu@INCResearch.com

•	Reference
HY Lin, DT Chen, PY Huang et al. SNP Interaction Pattern Identifier (SIPI): An Intensive Search for SNP-SNP Interaction Patterns. (under review).  

•	Example codes for Pairwise interaction analyses

******    Macro parameters    ******;
%let INdata = p.studpdata;   
%let OUTdt1 = p.AlleleMAF;   
%let OUTdt2 = p.BestInteract; 

%let SPvar1 = ;
%let SPvar2 = ;

%let OUTvar = d; ** binary outcome ;
%let Covariate_C = gender group;
%let Covariate_N = age;
%let OTHvar = ID; ** List of non-SNP, non-outcome variables;

******    Actions   ******;
%SIPI;

************/;


%macro SIPI ;
  proc format ;
    value LRLf 1 = '_Full     '
               2 = '_M1_int_o1'
               3 = '_M1_int_r1'
               4 = '_M2_int_o2'
               5 = '_M2_int_r2'
               6 = '_int_oo   '
               7 = '_int_ro   '
               8 = '_int_or   '
               9 = '_int_rr   ' ;
  run ;
  %global spv1 spv2 ;

  %macro LR9m (spv1, spv2);
     %do LRI = 1 %to 5 ;
        %if &LRI= 1 %then %let LRm= DD ;
        %else %if &LRI= 2 %then %let LRm= DR;
        %else %if &LRI= 3 %then %let LRm= RD;
        %else %if &LRI= 4 %then %let LRm= RR;
        %else %if &LRI= 5 %then %let LRm= AA;
        %do LRL = 1 %to 9 ;
           %if &LRL= 1 %then %let cov= %lowcase(%substr(&LRm,1,1))&spv1. | %lowcase(%substr(&LRm,2,1))&spv2. ;
           %else %if &LRL= 2 %then %let cov=  %lowcase(%substr(&LRm,1,1))&spv1.    %lowcase(%substr(&LRm,1,1))&spv1 *  %lowcase(%substr(&LRm,2,1))&spv2. ;
           %else %if &LRL= 3 %then %let cov= r%lowcase(%substr(&LRm,1,1))&spv1.   r%lowcase(%substr(&LRm,1,1))&spv1 *  %lowcase(%substr(&LRm,2,1))&spv2. ;
           %else %if &LRL= 4 %then %let cov=  %lowcase(%substr(&LRm,2,1))&spv2.    %lowcase(%substr(&LRm,1,1))&spv1 *  %lowcase(%substr(&LRm,2,1))&spv2. ;
           %else %if &LRL= 5 %then %let cov= r%lowcase(%substr(&LRm,2,1))&spv2.    %lowcase(%substr(&LRm,1,1))&spv1 * r%lowcase(%substr(&LRm,2,1))&spv2. ;
           %else %if &LRL= 6 %then %let cov=  %lowcase(%substr(&LRm,1,1))&spv1. *  %lowcase(%substr(&LRm,2,1))&spv2. ;
           %else %if &LRL= 7 %then %let cov= r%lowcase(%substr(&LRm,1,1))&spv1. *  %lowcase(%substr(&LRm,2,1))&spv2. ;
           %else %if &LRL= 8 %then %let cov=  %lowcase(%substr(&LRm,1,1))&spv1. * r%lowcase(%substr(&LRm,2,1))&spv2. ;
           %else %if &LRL= 9 %then %let cov= r%lowcase(%substr(&LRm,1,1))&spv1. * r%lowcase(%substr(&LRm,2,1))&spv2. ;

           data Fit ;  **  Create a dummy FIT dataset in case of no 'FitStatistics' output ;
               attrib Criterion length=$9 label='Criterion' format=$9.
                      InterceptAndCovariates format=9.3 ;
               Criterion = 'SC' ;
               InterceptAndCovariates = . ;
           run ;

           %put Model =&LRm., Seq #=&LRI.,  Cov #=&LRL.,  Cov= &Cov. ;
           proc logistic data=InheritanceMode ;
              %if %nrbquote(&covariate_c.) ne %then class &covariate_c. ;;
              model &OUTvar.(event='1' ref='0') = &cov. &covariate_c. &covariate_n. ;
              ods output ParameterEstimates=Param FitStatistics=Fit;
           run;

      ***  Collecting interaction chisq and p-value  *** ;
           data ParamCP (rename=(WaldChiSq=int_chisq ProbChiSq=int_p) drop=Variable);
              set Param (keep=Variable WaldChiSq ProbChiSq) ;
              where index(variable, '*')>0 ;
           run;

      ***  Collecting fits  *** ;
           data FitBIC (rename=(InterceptAndCovariates=BIC) drop=Criterion) ;
              set Fit (keep=Criterion InterceptAndCovariates) ;
              where Criterion='SC';
           run;

           data LR9m_single ;
              Length Var1 Var2 Model $40. ;
              if _n_ = 1 then set FitBIC ; set ParamCP ;
              Var1  = "&spv1" ;
              Var2  = "&spv2" ;
              Model = "&LRm." || put(&LRL, LRLF.) ;
              label BIC = "BIC" ;
           run;

           proc datasets library=work memtype=data nolist ;
              append base=Total data=LR9m_single force ;
              delete Param ParamCP fit FitBIC LR9m_single ;
           quit ;
        %end ;
     %end ;
  %mend ;

  %macro Action1 ;
     %if &spvp = 2 %then %do ;
         proc sort sortseq=uca (numeric_collation=on) data=Final out=&OUTdt1 ;
            by name ;
            where upcase(name) in ("%upcase(&SPvar1)" "%upcase(&SPvar2)") ;
      /* data _null_; set &OUTdt1;  call symput ('spv'||put(_n_, 1.), strip(name)); */
         run ;
         %let spv1 = &SPvar1 ;
         %let spv2 = &SPvar2 ;
         %PUT ===== One-pair SNP interaction analyses for %cmpres(&spv1) and %cmpres(&spv2) will be preformed ====== ;
         DM 'post "One-pair SNP interaction analyses for &SPvar1. and &SPvar2. will be preformed."';
         options notes ;
         ods listing ;
         ods results on ;
     %end ;
     %else %do ;
         proc sort sortseq=uca (numeric_collation=on) data=Final out=&OUTdt1 ;
            by name ;
         run ;
         %PUT ===== This program will run through all pair-wise SNP ====== ;
         %PUT ===== There will be %sysevalf(&spvs * (&spvs - 1) /2) combination of SNP pairs =====;
         DM 'post "Pairwise SNP interaction analyses will be performed."' ;
         options nonotes ;
         ods listing close ;
         ods html close ;
         ods results off ;
     %end ;
     proc datasets library=work memtype=data nolist ;
        delete allele xxfrq xxq_: varnms Major MAF Final ;
     quit ;

  %mend ;

  %macro Action2 ;
     %if &spvp = 2 %then %do ;
         %PUT ====== The program is running for just one SNP pair of %cmpres(&spv1) & %cmpres(&spv2) ====== ;
         %LR9m (&spv1, &spv2) ;
         data &OUTdt2._1u ; set Total ;
            if int_chisq = . then BIC = . ; ** BIC is set to null when no interactions ;
         proc sort data=&OUTdt2._1u out=&OUTdt2._1 ;
            by BIC int_chisq ;
         run ;
     %end ;
     %else %do ;
         %PUT ====== The program is running for all SNP pairs ====== ;
         proc transpose data=&OUTdt1. out=IntVars (drop=_name_);
            var name ;
         data IntCmbo ;
            set IntVars ;
            array SNPS (*) COL: ;
            do I = 1 to dim(SNPS)-1 ;
               do J = (I+1) to dim(SNPS) ;
                  SN1 = strip(SNPS (I)) ;
                  SN2 = strip(SNPS (J)) ;
                  output ;
               end ;
            end ;
         run ;
         data _null_ ;
            set IntCmbo end=eof ;
            call symput ('SPVa'||left(put(_n_, best.)), SN1) ;
            call symput ('SPVb'||left(put(_n_, best.)), SN2) ;
            if eof then call symput ('nCombo', left(put(_n_, best.))) ;
         run ;
         %do k = 1 %to &nCombo ;
            %PUT ====== The program is running on SNP pair of %cmpres(&&SPVa&k) & %cmpres(&&SPVb&k) ====== ;
            %LR9m (&&SPVa&k, &&SPVb&k) ;
         %end ;
        data Total ; set Total ;
            if int_chisq = . then BIC = . ; ** BIC is set to null when no interactions ;
         proc sort data=Total ;
            by Var1 Var2 BIC ;
            where BIC > . and int_chisq > . ;
         data Total ;
            set Total ;
            by Var1 Var2 BIC ;
            if first.Var2 ;
         proc sort sortseq=uca (numeric_collation=on) data=Total (drop=BIC) out=&OUTdt2._all ;
             by Var1 Var2 ;
         run ;
         proc datasets library=work memtype=data nolist ;
            delete IntVars IntCmbo ;
         quit ;
     %end ;
     proc datasets library=work memtype=data nolist ;
        delete InheritanceMode Total ;
     quit ;
     options notes ;
     ods listing ;
     ods html ;
     ods results on ;
  %mend ;

  ******    Initial Checking of the binary outcome variable    ****** ;
  proc sql noprint ;
     select count(*) into:BV0 from &INdata. where &OUTvar. = 0 ;
     select count(*) into:BV1 from &INdata. where &OUTvar. = 1 ;
     select count(*) into:BVa from &INdata. where &OUTvar. > . ;
  quit ;
  %let BV0 = &BV0 ;  %let BV1 = &BV1 ;  %let BVa = &BVa ;

  %if &BV0=0 or &BV1=0 or &BVa > %eval(&BV0 + &BV1) %then %do ;
      %PUT ===== The program stops because the outcome variable &OUTvar. does not meet the requirements ====== ;
      DM 'post "The program stops because the outcome variable &OUTvar. does not meet the requirements."';
  %end ;
  %else %do ;
     ods html close ; ods html ;
     ******    Getting SNP data    ****** ;
     proc contents data=&INdata. (drop=&OUTvar. &Covariate_C. &Covariate_N. &OTHvar.) out=varnms (keep=name type where=(type=2)) noprint ;
     data _null_ ;
        set varnms end=eof ;
        if _n_ = 1 then call execute ("proc freq data=&INdata.(where=(&OUTvar. in (0 1))) noprint;") ;
        call execute ("table "||strip(name)||" / missing out=xxq_"||strip(put(_n_, best.))||" (drop=percent);") ;
        if eof then call execute ("run;") ;
     data _null_ ;
        set varnms ;
        call execute ("data xxq_"||strip(put(_n_, best.))||" (keep=name value count); set xxq_"||strip(put(_n_, best.))||";") ;
        call execute ("length name $32. value $4.; name='"||strip(name)||"'; value="||strip(name)||"; run;") ;
     run ;

     data xxfrq ;
        set xxq_: ;
        where value ne '' ;
     proc sort data=xxfrq ;
        by name value ;
     run ;

     data xxfrq ;
        set xxfrq ;
        by name value ;
        retain C T A G 0;
        if first.name then do ;
           C=0; T=0; A=0; G=0;
        end ;
        do pose = 1 to 2 ;
           if substr(value, pose, 1)='C' then C = C + count;
           if substr(value, pose, 1)='T' then T = T + count;
           if substr(value, pose, 1)='A' then A = A + count;
           if substr(value, pose, 1)='G' then G = G + count;
        end ;
        if last.name ;
     run ;

     proc transpose data=xxfrq out=allele (rename=(_name_=allele)) prefix=vcnt;
        var C T A G;
        by name ;
     run ;

     proc sort data=allele ;
        by name descending vcnt1 ;
        where vcnt1 > 0 ;
     data allele ;
        set allele ;
        by name descending vcnt1 ;
        retain vcnt1r ;
        if first.name then do ;
           Major = 1 ;
           vcnt1r = vcnt1 ;
        end ;
        else do ;
           Major + 1 ;
           vcnt1r = vcnt1r + vcnt1 ;
           MAF = round((vcnt1 / vcnt1r), .01) ;
        end ;
        if Major < 3 ;
     proc transpose data=allele out=Major (drop=_name_ _label_) prefix=M;
        var allele ;
        by name ;
        id major ;
     run ;

     proc sort data=allele (keep=name MAF) out=MAF ;
       by name ;
       where MAF > . ;
     proc sort data=varnms (drop=type) ;
        by name ;
     data Final ;
        merge varnms (in=a) Major (rename=(M1=Major M2=Minor)) MAF ;
        by name ;
        if Major ne Minor and cmiss(Major, Minor)=0 ;
        Label Name  = "SNP Variables"
              Major = "Major Allele"
              Minor = "Minor Allele"
              MAF   = "Minor Allele Frequency" ;
     run ;

     ******     To determine 1-pair or all-pair outputs     ****** ;
     proc sql noprint ;
        select count(*) into:spvs
        from Final ;
        select count(*) into:spvp
        from Final (where=(upcase(name) in ("%upcase(&SPvar1)" "%upcase(&SPvar2)"))) ;
     quit ;
     %Action1 ;

     ******    Deriving SNP inheritance mode data based on major and minor allele     ****** ;
     data _null_;
        set &OUTdt1. end=eof;
        length A B $1.;
        A=Major;
        B=Minor;
        if _n_=1 then call execute ("data InheritanceMode (compress=yes); set &INdata.(where=(&OUTvar. in (0 1)));");
        call execute (     "if "||strip(name)||"='"||B||B||"' then r"||strip(name)||"=1;");
        call execute ("else if "||strip(name)||"='"||A||A||"'|"||strip(name)||"='"||A||B||"'|"||strip(name)||"='"||B||A||"' then r"||strip(name)||"=0;");
        call execute (     "if "||strip(name)||"='"||A||A||"' then d"||strip(name)||"=0;");
        call execute ("else if "||strip(name)||"='"||B||B||"'|"||strip(name)||"='"||A||B||"'|"||strip(name)||"='"||B||A||"' then d"||strip(name)||"=1;");
        call execute (     "if "||strip(name)||"='"||B||B||"' then a"||strip(name)||"=2;");
        call execute ("else if "||strip(name)||"='"||A||B||"'|"||strip(name)||"='"||B||A||"' then a"||strip(name)||"=1;");
        call execute ("else if "||strip(name)||"='"||A||A||"' then a"||strip(name)||"=0;");
        call execute (      "rr"||strip(name)||"=1-r"||strip(name)||";");
        call execute (      "rd"||strip(name)||"=1-d"||strip(name)||";");
        call execute (      "ra"||strip(name)||"=2-a"||strip(name)||";");
        call execute ("label  d"||strip(name)||"='"||strip(name)||" dominant b/on minor allele ' ");
        call execute (       "r"||strip(name)||"='"||strip(name)||" recessive b/on minor allele' ");
        call execute (       "a"||strip(name)||"='"||strip(name)||" additive b/on minor allele ' ");
        call execute (      "rd"||strip(name)||"='"||strip(name)||" dominant b/on major allele ' ");
        call execute (      "rr"||strip(name)||"='"||strip(name)||" recessive b/on major allele' ");
        call execute (      "ra"||strip(name)||"='"||strip(name)||" additive b/on major allele ';");
        if eof then call execute ("run;");
     run ;
     %Action2 ;
  %end ;
%mend ;
