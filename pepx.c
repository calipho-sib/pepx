// For now update both at git github.com/calipho-sib/pepx and svn http://trac/svn/NextProt/main/trunk/com.genebio.nextprot.tools.sequence
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
//#include <sys/syscall.h>
#include <errno.h>
#include <getopt.h>
#include <malloc.h>

#define FALSE 0
#define TRUE  1
#define MAXSEQ 32768
// MAXSEQ max number of results for a query
#define MAXHEADSIZE 1500000
// MAXHEADSIZE size of titin longest isoform header line in PEFF
#define MAXSEQSIZE 36500
// MAXSEQSIZE = size of titin longest isoform
#define MAXISO 100
#define MAXVARINPEP 127
#define MAXVARPERAA 8
#define LINELEN 1024
#define ACLEN 16
//#define BINSIZE 32  // adjust w stats
#define BINSIZE 16  // adjust w stats, this is becoming a problem, have to reduce to 12 when indexing with variants !
// optimal binsizes: 2->1000, 3->100, 4->25, 5->10, 6->?
#define MINPEPSIZE 3
//#define MAXPEPSIZE 5
#define MAXPEPSIZE 6
// Allowing 7 implies adapting BINSIZE t0 low value in 7-mers (memory pbs)
#define SILENT 1
#define INDEX2WIDTH 12
#ifndef max
	#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

// Globals

typedef struct idxdata
       {
       FILE*  fh;
       FILE*  ffh;
       int elemcnt;
       } IDXDATA;

IDXDATA idxinfo[MAXPEPSIZE], idxinfoIL[MAXPEPSIZE];
FILE* currentindex;
char currentresult[32];
static int currentpepsize; // weird but alters idxinfo[6] when changed if not static

typedef struct bin
       {
       int  ac_cnt;
       struct bin* nextbinptr;
       char ac[BINSIZE][ACLEN];
       } BIN;

BIN** toplevel[MAXPEPSIZE+1];

int varindex_overflow_cnt = 0;
char variants[MAXSEQSIZE][MAXVARPERAA];
char results[MAXSEQ][ACLEN];
int aacode[26] = {0,2,1,2,3,4,5,6,7,7,8,9,10,11,-1,12,13,14,15,16,1,17,18,-1,19,3};
// B is D, Z is E, U is C, J is I
char aarevcode[] = "ACDEFGHIKLMNPQRSTVWY";
char aarevcodeIL[] = "ACDEFGHJKLMNPQRSTVWY";
char matchmode[10] = "ACISO";
int pow20[] = {1,20,400,8000,160000,3200000,64000000,1280000000};
int binsizes[] = {0,0,1000,100,25,10,6,4};
char masterseq[MAXSEQSIZE];
char currAC[ACLEN];
char currISO[ACLEN];
char outputmode[16];
char linesep[4] = "\n";
char envstring[LINELEN];
char version[] = "1.57";
// debug/profiling/stats stuff
int debug;
int totalbins=0;
clock_t start;
double duration;
int umatchcnt;
int maxoccur=0;
int usscnt= 0;
int totalvarcnt = 0;
int varmisscnt1 = 0, varmisscnt2 = 0, varmisscnt3=0;
int ignore_variants = 0;
int IL_merge = 0;
int json = 0;
int FASTAinput = 0;
static char jsonBuffer[8192];// weird but alters idxinfo[6] when changed if not static
char varfolder[LINELEN];
char indexfolder[LINELEN]=".";

// ---------------- code4tenAA ---------------------
char* code4tenAA(char* currentAC)
{
int i, found=0;
static char codedAC[16]="";
//static char* tab[]={"A A0A087","B A0A0B4","C A0A075","D A0A096","E A0A0C5", "F A0A0J9","G A0A0C4", "H A0A0A0", "I A0A0G2"};
static char* tab[]={"A A0A024","B A0A044","C A0A067","D A0A075","E A0A087","F A0A088","G A0A096","H A0A0A0","I A0A0A6","J A0A0B4","K A0A0C3",
          "L A0A0C4","M A0A0C5","N A0A0D2","O A0A0D5","P A0A0D9","Q A0A0E3","R A0A0G2","S A0A0H2","T A0A0J9","U A0A0K0","V A0A0M3","W A0A0N4",
          "X A0A0R4","Y A0A0S2","Z A0A0U1","a A0A0U5","b A0A0X1","c A0A140","d A0A182","e A0A1B0","f A0A1C7","g A0A1D5"};

if(!strncmp(currentAC,"P",1))
  // convert back to 10 digit AC
  {
  strncpy(codedAC,currentAC+1,1);
  codedAC[1] = 0;
  for(i=0;i<sizeof(tab)/8;i++)
     if(!strncmp(tab[i],codedAC,1)) 
      {
      found=1;
      break;
      }
  if(found)
    {
    strcpy(codedAC,tab[i]+2);
    strcat(codedAC,currentAC+2);
    //fprintf(stderr,"output: %s\n",codedAC);
    return(codedAC);
    }
  else
    strcpy(codedAC,"XXXXXX");
  }
else
  // convert to 6 digit AC
  {
  strncpy(codedAC,currentAC,6);
  codedAC[6]=0;
  for(i=0;i<sizeof(tab)/8;i++)
    if(strstr(tab[i],codedAC))
      {
      found=1;
      break;
      }
   if(found)
    {
    strcpy(codedAC,"P");
    strncat(codedAC,tab[i],1);
    strcat(codedAC,currentAC+6);
    }
  else
    strcpy(codedAC,"XXXXXX");
  }    
//fprintf(stderr,"output: %s\n",codedAC);
return(codedAC);
}
// ---------------- code2sseq ---------------------
void code2sseq(int sseqcode, int pepsize, char* subseq)
{
char sseq[16]="";
int i, j;

for(i=pepsize-1;i>=0;i--)
   {
   j = sseqcode / pow20[i];
   sseqcode %= pow20[i];
   if(IL_merge)
     strncat(sseq,&aarevcodeIL[j],1);
   else
     strncat(sseq,&aarevcode[j],1);
   }
strcpy(subseq,sseq);
}

// ---------------- sscode ---------------------
int sscode(char* subseq)
{
int sseqcode = 0, len, i;
len = strlen(subseq);
for(i=len;i>0;i--)
  sseqcode += pow20[len-i] * (aacode[subseq[i-1] - 'A']);
return(sseqcode);
}

// ---------------- pepx_json_header---------------------
void pepx_json_header(char* title, char* querystring, int modeIL)
// recapitulates given parameters
{

sprintf(jsonBuffer,"\"%s\":{\n",title);
fputs(jsonBuffer,stdout); // outputs header title
if(querystring != NULL)
  {
  sprintf(jsonBuffer,"\"peplist\":\"%s\",\n\"modeIL\": %d\n},",querystring,modeIL);
  fprintf(stdout,"%s",jsonBuffer);
  //fprintf(stderr,"jsonBuffer: %s\n",jsonBuffer);
  }
}

// ---------------- pepx_tojson ---------------------
char* pepx_tojson(char* finalresult)
{
int varpos=0;
char* dashptr;

//fprintf(stdout,"pepx_tojson receives %s\n",finalresult);
dashptr=strrchr(finalresult,'-');
if((dashptr != finalresult + 6) &&  (dashptr != finalresult + 10)) // don't forget 10 digirs
  // we have a variant pos in the last digits (eg: NX_Q8IWA4-2-121)
  {
  varpos = atoi(dashptr + 1);
  *dashptr = 0; // remove parsed pos from string
  }
sprintf(jsonBuffer,"{\n\"isoformName\":\"%s\",\n\"position\":%d\n}",finalresult, varpos);
return(jsonBuffer);
}

// ---------------- pepx_create_variant_files ---------------------
int pepx_create_variant_files(char* varfolder)
{
FILE *varmaster, *out;
char *token, buf[LINELEN], entryvarbuf[524288]="", org[64], dest[128], fname[128], lastac[16]="", ac[16], nxid[16], annotId[32];
int startpos, endpos, i=0;

fprintf(stderr,"Building for first time with variants, please wait...\n");
sprintf(fname,"pepx-variantdata.csv");
if((varmaster=fopen(fname,"r"))==NULL)
 {
 perror(fname);
 return(0);
 }
while(fgets(buf,LINELEN,varmaster))
     // Sample line: "NX_P16333-1,AN_P16333_000975,183,183,V,I", actually commas are replaced by spaces
     {
     if(strncmp(buf,"NX_",3))
       continue;
     *dest=0; //To avoid that last value is keeped when scaning a miss variant (dest empty)
     sscanf(buf,"%s %s %d %d %s %s\n",nxid, annotId, &startpos,&endpos,org,dest);
     //fprintf(stderr,"read %s\n",buf);
     //fprintf(stderr,"nxid: %s\n",nxid);
     strcpy(ac,nxid+3);
     if(strlen(lastac) && strcmp(ac,lastac))
       // write previous iso record 
       {
       sprintf(fname,"%s/%s.csv",varfolder,lastac); 
       //fprintf(stderr,"writing %s\n",fname);
       out=fopen(fname,"w");
       fprintf(out,"%s",entryvarbuf);
       fclose(out);
       if((++i % 2500) == 0)
         fprintf(stderr,"%d isoform variant files created...\n",i);
       sprintf(buf,"%d %d %s %s\n",startpos,endpos,org,dest);
       strcpy(entryvarbuf,buf);
       strcpy(lastac,ac);
       }
     else
       {
       sprintf(buf,"%d %d %s %s\n",startpos,endpos,org,dest);
       strcat(entryvarbuf,buf);
       strcpy(lastac,ac);
       }
     }
// write last record
sprintf(fname,"%s/%s.csv",varfolder,ac); 
out=fopen(fname,"w");
fprintf(out,"%s",entryvarbuf);
fclose(out);

sprintf(fname,"%s/VARIANTS_OK",varfolder);
out=fopen(fname,"w");
fclose(out);
return(1);
}

// ---------------- pepx_build_varindex ---------------------
int pepx_build_varindex(char* isoname)
{
char buf[LINELEN], orgAA[16], lastorgAA[16], varAA[16], fname[64], *ptr;
char iso[16];
int i, skipvar = 0, cmdres, fpos, lastfpos, lpos, localvarcnt=0;
FILE *NextProtvariants;

if(debug)
  fprintf(stderr,"Building variants for %s\n",currISO);
if(is10digits(isoname)) 
  // get around 10-len accs
  {
  strcpy(iso,code4tenAA(isoname));
  }
else strcpy(iso,isoname);

// Get variants from NextProt
sprintf(fname,"%s/%s.csv",varfolder,iso);
if(debug) fprintf(stderr,"opening file: %s\n",fname);
if((NextProtvariants=fopen(fname,"r"))==NULL)
  // THe variant file doesn't exist for this iso : no variants
  return(0);
  
// Reset table each time  
memset(variants,0,MAXSEQSIZE*MAXVARPERAA);
// Todo: check that orgAA is exactly the same (seq/vars desynchronization) 
while(fgets(buf,LINELEN,NextProtvariants))
     {
     *varAA=0; //To avoid that last value is keeped when scaning a miss variant (varAA empty)
     if(debug) fprintf(stderr,"variants line: %s\n",buf);
     sscanf(buf,"%d %d %s %s\n",&fpos,&lpos,orgAA, varAA);
     //if(strstr(iso,"P04637") && strstr(buf,"description") && strstr(buf,"sporadic"))
     if(fpos==1 || strchr(varAA,'*') || (strlen(orgAA) > 2)  || (strlen(varAA) > 2))
        // No variants on Met-1 allowed // No stop variants allowed
       {
       if((strlen(varAA) == 1 && orgAA[strlen(varAA)-1] == varAA[0]) || strlen(varAA) == 0)
	     varmisscnt3++;
       skipvar = 1;
       }
     if(skipvar)
       {
       skipvar = 0;
       continue;
       }
     if((fpos != lpos)) // not a simple snp
       {
       if((strlen(orgAA) == 2) && (strlen(varAA) == 2) && orgAA[0] == varAA[0])// like A7XYQ1(180-181)
	   // Store variant at 2nd pos
	 {
   if(IL_merge && (varAA[1] == 'I' || varAA[1] == 'L') )
    varAA[1] = 'J';
	 variants[fpos][0]=varAA[1];
	 localvarcnt++;
	 }
       else if((strlen(orgAA) == 2) && (strlen(varAA) == 1) && varAA[0] == orgAA[1]) // like A6NFR6-3
	   // Store variant
	   {
	   if((fpos == lastfpos+1) && (orgAA[0] == lastorgAA[0])) // To avoid damage like in O00555-1 2219-2220 and 2220-2221
	     {
	     lastfpos = fpos;
	     continue;
	     }
	   variants[fpos-1][0]='X'; // Missing
	   varmisscnt2++;
	   strcpy(lastorgAA,orgAA);
	   lastfpos = fpos;
	   localvarcnt++;
	   }
	 } 
       else 
	 {
	 if(IL_merge && ((orgAA[0] == 'I' && varAA[0] == 'L') || (orgAA[0] == 'L' && varAA[0] == 'I'))) // In IL mode variants I->L or L->I would just add noise and are ignored
	   continue;
	 // Store variant
	 if(strlen(varAA) == 0)
	   // code for 'Missing' AA
	   {
	   strcpy(varAA,"X");
	   varmisscnt1++;
	   }
	 for(i=0;i<MAXVARPERAA;i++)
	    // we allow 8 different variants at each position: Check if we already have a variant at this pos
	    if(variants[fpos-1][i] == 0)
	      break;
	 //if(i >= MAXVARPERAA) fprintf(stderr,"%s: multiple variants at %d-%d\n",currISO,fpos,lpos); 
   if(IL_merge && (varAA[0] == 'I' || varAA[0] == 'L') )
    varAA[0] = 'J';
	 variants[fpos-1][i]=varAA[0];
	 localvarcnt++;
	 }
       }

fclose(NextProtvariants);
if(debug)
  fprintf(stderr,"Variants done.\n");
  
return(localvarcnt);
}

// ---------------- pepx_loadall ---------------------
void pepx_loadall()
{
int i, pepcnt, filelen;
char fname[LINELEN], *ptr, idxpath[128]="";
FILE* idx;

if(strlen(indexfolder) > 8)
  {
  strcpy(idxpath,indexfolder);
  strcat(idxpath,"/");
  }
else if((ptr=getenv("PEPX")) != NULL)
  {
  strcpy(idxpath,ptr);
  strcat(idxpath,"/");
  }
for(i=MINPEPSIZE; i<= MAXPEPSIZE;i++)
   {
   usscnt = 0;
   if(IL_merge)
     sprintf(fname,"%spepxIL%d.idx",idxpath,i);
   else  
     sprintf(fname,"%spepx%d.idx",idxpath,i);
    if((idx=fopen(fname,"r"))==NULL)
     {
     perror(fname);
     exit(2);
     }
   // save flat file handle
   if(IL_merge)
      idxinfoIL[i].ffh = idx;
   else
      idxinfo[i].ffh = idx;
   strcat(fname,"2");
   if((idx=fopen(fname,"r"))==NULL)
     {
     perror(fname);
     exit(2);
     }
   fseek(idx,0,SEEK_END);
   filelen = ftell(idx);
   pepcnt = filelen/(i+INDEX2WIDTH);     
   //printf("%d-mers size: %d bytes -> %d elts\n",i,filelen,pepcnt);
   if(IL_merge)
     {
     idxinfoIL[i].fh = idx;
     idxinfoIL[i].elemcnt = pepcnt;
     }
   else
     {
     idxinfo[i].fh = idx;
     idxinfo[i].elemcnt = pepcnt;
     }
   }
}

// ---------------- pepx_initindexes ---------------------
void pepx_initindexes()
{
int pepsize;
BIN **ptr;

for(pepsize=MINPEPSIZE; pepsize <= MAXPEPSIZE;pepsize++)
   {
   if((ptr = calloc(1, pow20[pepsize] * sizeof(BIN*))) == NULL)
     {
     fprintf(stderr,"Allocation failed\n");
     exit(11);
     }

   //printf("%d empty ints reserved at address %p\n",pow20[pepsize],ptr);
   toplevel[pepsize] = ptr;
   }
}

// ---------------- pepx_save ---------------------
void pepx_save(char* fname, int pepsize)
{
int sseqcode, i, maxcode, ac_cnt, iso_cnt;
char ac[ACLEN], ac_noiso[ACLEN], lastac[ACLEN], subseq[16], maxseq[16]="";
FILE *out, *out2;
BIN *binptr, **binarray;
int idxpos;

if((out=fopen(fname,"w"))==NULL)
 {
 perror(fname);
 exit(2);
 }

// secundary index
strcat(fname,"2");
if((out2=fopen(fname,"w"))==NULL)
 {
 perror(fname);
 exit(2);
 }
// index header line (1rst pep line must not be at pos 0 for new index schema)
for(i=0;i<pepsize+INDEX2WIDTH-1;i++)
   fprintf(out2,"0");
fprintf(out2,"\n");

usscnt = maxoccur = umatchcnt = 0;
fprintf(stderr,"Writing index for %d-mers ...\n",pepsize);
binarray = toplevel[pepsize];
maxcode = pow20[pepsize];
for(sseqcode=0; sseqcode<maxcode; sseqcode++)
   if((binptr=binarray[sseqcode]) != NULL)
     // This pep has some matches
     {
     code2sseq(sseqcode, pepsize, subseq);
     idxpos = ftell(out);
     if(fprintf(out,"%s\n",subseq) <= 0)
       fprintf(stderr,"writing acs for %s failed\n",subseq);
     // write secundary index
     fprintf(out2,"%s %10u\n",subseq, idxpos);
     ac_cnt = iso_cnt = 0;
     usscnt++;
     while(binptr->nextbinptr != NULL)
          {
	  // add the next BINSIZE ac
	  for(i=0;i<BINSIZE;i++)
	     {
	     strcpy(ac,binptr->ac[i]);
	     fprintf(out,"%s\n",ac);
	     //fprintf(stderr,"ac: %s\n",ac);
	     iso_cnt++;
	     strcpy(ac_noiso,ac);
	     ac_noiso[6] = 0;
	     if(ac_noiso != lastac)
	       ac_cnt++;
	     strcpy(lastac,ac_noiso); 
	     }
	  // set pointer to next bin
	  binptr = binptr->nextbinptr;
	  }
     // dump last bin
     for(i=0;i<binptr->ac_cnt;i++)
	{
	strcpy(ac,binptr->ac[i]);
	fprintf(out,"%s\n",ac);
        iso_cnt++;
        strcpy(ac_noiso,ac);
        ac_noiso[6] = 0;
        if(ac_noiso != lastac)
	  ac_cnt++;
	strcpy(lastac,ac_noiso);
	}
	
     if(ac_cnt == 1)
       umatchcnt++;
     else if(ac_cnt > maxoccur)
       {
       maxoccur =  ac_cnt;
       strcpy(maxseq, subseq);
       }			    
     }
fprintf(stderr,"%d distincts %d-mer, %d are uniques (found in just 1 entry), most frequent: %s (%d)\n",usscnt, pepsize, umatchcnt, maxseq, maxoccur);
fclose(out);
fclose(out2);
}

// ---------------- pepx_saveall ---------------------
void pepx_saveall()
{
int i;
char fname[LINELEN];

for(i=MINPEPSIZE; i<= MAXPEPSIZE;i++)
   {
   if(!IL_merge)  
     sprintf(fname,"%s/pepx%d.idx",indexfolder,i);
   else
     sprintf(fname,"%s/pepxIL%d.idx",indexfolder,i);
   pepx_save(fname, i);
   }
}

/******* rescompare **************************************************/  

int rescompare (char *str1, char *str2)    
{
return (strncmp(str1,str2,ACLEN));
}

/******* pepcompare **************************************************/  

int pepcompare (char *str1, char *str2)    
{ 
fseek(currentindex,(long)str2,0);
fgets(currentresult,32,currentindex);
//fprintf(stderr,"current idx2 file pos:  %li, fseek returned %d\n",(int)str2,fret);
//fprintf(stderr,"res: %s\n",currentresult);
return (strncmp(str1,currentresult,currentpepsize));
}

// ---------------- pepx_reportnomatch ---------------------
int pepx_reportnomatch(char* orgquerystring)
{
//fprintf(stdout,"NO_MATCH %s\n", outputmode);
if(json)
  fprintf(stdout,"]\n}"); // close empty entrymatches
else  
  fprintf(stdout,"NO_MATCH %s%s", orgquerystring,linesep);
return(0);
}

// ---------------- pepx_merge_with_prev_res ---------------------
int pepx_merge_with_prev_res(char endres[MAXSEQ][ACLEN+4], char curres[MAXSEQ][ACLEN], char* acstr, int rescnt)
{
char isoonly[ACLEN], isovar[ACLEN]="", prevmatch[ACLEN]="", *dashptr, *matchptr;  
int i, j, prevmatchvar, dupcnt=0;

//fprintf(stderr,"acstring: %s\n",acstr);
j = rescnt;
if(strlen(acstr) == 0)
  // first pep of a query
  {
  //fprintf(stderr,"\nfirst 6-mer: %d hits\n",j);  
  for(i=0;i<rescnt;i++)
     {
     if(i && !strcmp(curres[i],curres[i-1]))
     {
     dupcnt++;
     continue; // fprintf(stderr,"\ndup: %s\n",curres[i]);, TODO: investigate why it occurs, eg: ASQSIS
     }
     //fprintf(stderr,"%d: %s (dupcnt: %d)\n",i,curres[i],dupcnt);
     strcpy(endres[i],curres[i]);
     strcat(acstr,curres[i]);
     strcat(acstr,",");
     }
   j = rescnt - dupcnt;  
  }   
else // Check the match did exist for previous peps
  {
  for(i=j=0;i<rescnt;i++)
     {
     strcpy(isoonly,curres[i]);
     if((dashptr=strrchr(isoonly,'-')) != isoonly + 6)
       // Remove variant pos before searching prev results
       {
       //fprintf(stderr,"removing varpos for %s\n",isoonly);
       strcpy(isovar,isoonly);
       //fprintf(stderr,"isovar now: %s\n",isovar);
       *dashptr = 0;
       }
     //fprintf(stderr,"checking %s (result %d/%d)\n",curres[i],i+1,rescnt);
     if(matchptr=strstr(acstr,curres[i]))
       for(;matchptr;matchptr=strstr(matchptr+1,curres[i]))
       {
       strncpy(prevmatch,matchptr,ACLEN);
       *strchr(prevmatch,',')=0;
       if((dashptr=strrchr(prevmatch,'-')) != prevmatch + 6)
	 prevmatchvar = TRUE;
       else
	 prevmatchvar = FALSE;
	 //fprintf(stderr,"prevmatchvar: %d\n",prevmatchvar);
       if(!strcmp(prevmatch,curres[i]) || (prevmatchvar && strstr(prevmatch,curres[i])))
	 // If a new match without variant falls on a previous match with variant: keep variant
         {
         //fprintf(stderr,"keeping prevmatch: %s\n",prevmatch);
         strcpy(endres[j++],prevmatch);
	 }
       //else fprintf(stderr,"discarding prevmatch: %s\n",prevmatch);	 
       }
     else if (isovar &&  (matchptr=strstr(acstr,isoonly)))
       // new match is a variant
       {
       strncpy(prevmatch,matchptr,ACLEN);
       *strchr(prevmatch,',')=0;
       //fprintf(stderr,"isovar: %s, prevmatch: %s\n",isovar,prevmatch);
       if(strrchr(prevmatch,'-') == prevmatch + 6)
         // replace simple match by new variant match
         strcpy(endres[j++],isovar);
       else if(strstr(prevmatch,isoonly))
	     // reconduct previous variant pos
      strcpy(endres[j++],prevmatch);
       //else fprintf(stderr,"%s not added to endres\n\n",isovar);
       }
     *isovar = 0;  
     }
     
  if(j)  
   // regenerate actring
   {
   //fprintf(stderr,"j equals %d, regenerating ACstring\n",j);
   *acstr = 0;
   for(i=0;i<j;i++)
      {
      strcat(acstr,endres[i]);
      strcat(acstr,",");
      }
    //fprintf(stderr,"regenerated ACstring: %s\n\n",acstr);  
    }
  }
return(j);
}

// ---------------- pepx_search ---------------------
int pepx_search(char* query, IDXDATA* idx)
{
typedef char actab[ACLEN];  
int i, j, jpos=-1, jcnt=0, found, rescnt=0, pepsize, fpos;
char querystring[MAXPEPSIZE], ac[ACLEN], newquery[MAXPEPSIZE], fname[LINELEN], buf[MAXPEPSIZE+1], *ptr;
char acholder[ACLEN+1];
actab *currentresults;
FILE* ffh;

strcpy(querystring,query);
pepsize = strlen(query);
currentpepsize = pepsize;
currentindex = idx[pepsize].fh;
currentresults = results;

for(i=0;i<pepsize;i++)
   if(querystring[i] == 'X')
     {
     jpos = i;
     jcnt++;
     }
   else if(IL_merge && (querystring[i] == 'I' || querystring[i] == 'L'))
     // We are querying against the IL_merged index
     querystring[i] = 'J';
     
if(jcnt > 1)
  // joker
  {
  fprintf(stderr,"\n%s: No more than 1 joker/6AA\n",query);
  return(0);
  }
else if((jcnt == 1) && (jpos != 0) && (jpos != pepsize-1))
  // internal joker
  {
  for(i=0;i<strlen(aarevcode);i++)
     // For each of the 20 AAs
     {
     memset(newquery,0,MAXPEPSIZE);
     strncpy(newquery,querystring,jpos);
     strncat(newquery,&aarevcode[i],1); // Todo: check if this works for I/L
     strcat(newquery,querystring + jpos + 1);
     //fprintf(stderr,"\nnewquery: %s\n",newquery);
     if(bsearch(newquery,NULL,idx[pepsize].elemcnt,pepsize+INDEX2WIDTH,pepcompare))
       {
       fpos = atoi(currentresult + pepsize);
       ffh = idx[pepsize].ffh;
       fseek(ffh,fpos,SEEK_SET);
       // skip peptide
       fgets(ac,ACLEN,ffh);
       while(fgets(ac,ACLEN,ffh))
            // get all ACs
            {
            if(strlen(ac) > 7)
              {
              if(ptr=strrchr(ac,'\n'))
	        *ptr = 0;
              else // max len identifiere eg:Q8WZ42-11-11658
                ac[ACLEN-1] = 0;
              // make results uniques
	      for(found=FALSE,j=0;j<rescnt;j++)
	         if(!strcmp(ac,currentresults[j]))
		   {
		   found = TRUE;
		   break;
		   }
	      if(!found)
                strcpy(currentresults[rescnt++],ac);
              }
            else
              break;
            }
       } 
     }
  return(rescnt);
  }

if(jpos == 0)
  // initial  joker -> reduce query length
  {
  strcpy(buf,querystring+1);
  strcpy(querystring,buf);
  currentpepsize = --pepsize;
  currentindex = idx[pepsize].fh;
  }
else if (jpos == pepsize-1)
  // terminal  joker -> reduce query length
  {
  querystring[pepsize-1] = 0;
  currentpepsize = --pepsize;
  currentindex = idx[pepsize].fh;
  }
  
// No jokers
if(!bsearch(querystring,NULL,idx[pepsize].elemcnt,pepsize+INDEX2WIDTH,pepcompare))
  // No match
  {
  //fprintf(stderr,"\n%s: not found in %s\n",querystring,idx==idxinfoIL?"ILindex":"BaseIndex");  
  return(0);
  }
//fprintf(stderr,"\n%s: found in %s\n",querystring,idx==idxinfoIL?"ILindex":"BaseIndex");
fpos = atoi(currentresult + pepsize);
ffh = idx[pepsize].ffh;
fseek(ffh,fpos,SEEK_SET);
// skip peptide
fgets(ac,ACLEN,ffh);
//fprintf(stderr,"\nReading match ACs for %s\n",ac);
while(fgets(acholder,ACLEN+1,ffh))
     // get all ACs
     {
     strncpy(ac,acholder,ACLEN);  
     if(strlen(ac) > 7)
       {
       if(ptr=strrchr(ac,'\n'))
	 *ptr = 0;
       else // max len identifiere eg:Q8WZ42-11-11658
         ac[ACLEN-1] = 0;
       //fprintf(stderr,"\nrescnt %d %s: %s(%d)\n",rescnt, query,ac,strlen(ac));
       strcpy(currentresults[rescnt++],ac);
       }
     else
       {
       //fprintf(stderr,"\nBreaking at %s\n",ac);	 
       break;
       }
     }
//fprintf(stderr,"\nReturning %d matches\n",rescnt);     
return(rescnt);
}

int is10digits(char* match)
{
if(match[0] == 'P' && isalpha(match[1]))  
/*if(!strncmp(match,"PA",2) || !strncmp(match,"PB",2) || !strncmp(match,"PC",2) ||
   !strncmp(match,"PD",2) || !strncmp(match,"PE",2) || !strncmp(match,"PF",2) ||
   !strncmp(match,"PG",2) || !strncmp(match,"PH",2) || !strncmp(match,"PI",2)) */
 return(1);
return(0);  
}

// ---------------- pepx_displaymatches ---------------------
void  pepx_displaymatches(char matches[][ACLEN+4] , int nbmatch, char* orgquerystring, char* querystring)
// orgquerystring the original modified or commented peptide (GCS*PLKKT*V), querystring the peptide stripped on non AA characters
{
int ac_cnt, i;
char ac[ACLEN], lastac[ACLEN]="";
  
// Sort results by AC
if(nbmatch > 2)
  qsort(matches, nbmatch, ACLEN+4, rescompare);
  
if(!strcmp(matchmode,"ACONLY"))
   // Re-filter without iso ids (it is already sorted)
   {
   ac_cnt = 0;
   for(i=0;i<nbmatch;i++)
      {
      if(!strchr(matches[i],'-')) 
        // it can happen at pos 0, eg: VIDLT, to investigate
        {
        //fprintf(stderr,"No dash at %d %s\n",i,ac);
        continue;
        }
      strcpy(ac,matches[i]);
      *strchr(ac,'-') = 0;
      if(strcmp(ac,lastac)) 
	       strcpy(results[ac_cnt++],ac);
      strcpy(lastac,ac);
      }
   // re-fill final results
   nbmatch = ac_cnt;   
   for(i=0;i<nbmatch;i++)
      strcpy(matches[i],results[i]);
   }  
   
if(!strcmp(outputmode,"BATCH"))
  // input from file (results are comma-separated)
  {
  //printf(stdout,"\n");
  if(nbmatch == 0)
    fprintf(stdout,"NO_MATCH");
  else for(i=0;i<nbmatch;i++)
     {
     // Check for coded 10-digit ACs
     fprintf(stdout,"%s",is10digits(matches[i])?code4tenAA(matches[i]):matches[i]);
     if(i != nbmatch-1)
       fprintf(stdout,","); // match separator
     }
  // Reminder of the query
  fprintf(stdout," %s\n",orgquerystring);
  }
else
  {
  // input from stdin or envstring in web cgi mode (results are crlf-separated)
  //fprintf(stdout,"\n%s: %d match(s)%s",querystring,cnt,linesep);
  if(!json)
    fprintf(stdout,"%s%s: %d match(s)%s",linesep,querystring,nbmatch,linesep);
  for(i=0;i<nbmatch;i++) 
     {
       if(json)
         {
	 strcpy(ac,matches[i]);
         *strchr(ac,'-')=0;
	 if(strcmp(ac,lastac))
	   {
	   if(i)
	     // close previous iso and entrymatch tag
	     fprintf(stdout,"]\n},");
	   if(is10digits(ac)) // Check for coded 10-digit ACs
	     fprintf(stdout,"{\n\"entryName\":\"%s\",\n\"isoforms\":[",code4tenAA(ac));
	   else  
	     fprintf(stdout,"{\n\"entryName\":\"%s\",\n\"isoforms\":[",ac);
	   strcpy(lastac,ac);
	   }
	 else fprintf(stdout,",");
	 if(is10digits(ac))
	   fprintf(stdout,"%s",pepx_tojson(code4tenAA(matches[i])));
	 else
	   fprintf(stdout,"%s",pepx_tojson(matches[i]));
         }
       else fprintf(stdout,"%s%s",is10digits(matches[i])?code4tenAA(matches[i]):matches[i],linesep);
     }
    if(json)
      fprintf(stdout,"]\n}\n]\n}\n"); // finished, close iso and entrymatch tags tags 
    else
      fprintf(stdout,"%s: %d match(s)%s",querystring,nbmatch,linesep);
  }
}

// ---------------- pepx_filterquery ---------------------
int pepx_filterquery(char* query, char* querystring)
{
char *qptr;
int i;

for(i=0,qptr=querystring;i<strlen(query);i++)
   {
   // uppercase
   query[i] = toupper(query[i]);
   if(!strchr(aarevcode, query[i]))
     {
     if(isspace(query[i]) || isdigit(query[i])
     || (query[i] == '*')
     || (query[i] == '#')
     || (query[i] == '.')
     || (query[i] == '(')
     || (query[i] == ')'))
       // just don't copy unwanted, and allow to keep eventually tagged seq for result display
       i = i;
     else if((query[i] == 'X') || (query[i] == 'J' && IL_merge))
       // This is OK
       *qptr++ = query[i];
     else
       {
	 // Todo: Maybe consider U as C
       fprintf(stderr,"\n%c is not an AA\n",query[i]);
       return(0);
       }
     }
   else
     *qptr++ = query[i];
   }
*qptr=0;
return(1);
}

// ---------------- pepx_processquery ---------------------
int pepx_processquery(char* orgquerystring)
{
char query[8192], querystring[8192], nextquery[8192], subquery[MAXPEPSIZE + 1]="", acstring[400000]="", finalres[MAXSEQ][ACLEN+4];
char *qptr;
int row=0, i, j, cnt;

// TODO: rewrite and iterate one by one AA, or 1rst and last and one by one ?
strcpy(query,orgquerystring);
if(json) fprintf(stdout,"{\n\"peptide\":\"%s\",\n\"entryMatches\":[",query);
if(!strcmp(outputmode,"BATCH"))
  // Only first token is the peptide
  {
  for(i=0;i<strlen(query);i++)
     if(isspace(query[i]))
       {
       query[i] = 0;
       break;
       }
  }

// quick filter
 if(!pepx_filterquery(query, querystring))
   // invalid query
   return(0);
strcpy(query,querystring);
//if(strlen(query) >= 4096) {strncpy(subquery,query,12); subquery[12] = 0; fprintf(stderr,"\nQuery length: %d %s...\n",strlen(query),subquery);}
  
while(strlen(query) > MAXPEPSIZE)
   { // Split query in overlaping subqueries of maximum length 
   //fprintf(stderr,"pepx_processquery filtered query: %s, copying %d AAs in subquery\n",query,MAXPEPSIZE);
   strncpy(subquery,query,MAXPEPSIZE);
   subquery[MAXPEPSIZE] = 0;
   //fprintf(stderr,"pepx_processquery subquery %s\n",subquery);
   // Prepare next subquery
   strcpy(nextquery,query + MAXPEPSIZE-1); // Faster but more false pos
   //if(strlen(querystring) > 25)
     //strcpy(nextquery,query + MAXPEPSIZE-1); // Faster but more false pos
   /*else if(strlen(querystring) > 16)
     strcpy(nextquery,query + 2);
   else
     strcpy(nextquery,query + 1); */
   strcpy(query, nextquery);
   //fprintf(stderr,"will look for %s\n",query);
   if(IL_merge)
     cnt = pepx_search(subquery,idxinfoIL);
   else
     cnt = pepx_search(subquery,idxinfo);
   if(cnt==0)
     {
     //fprintf(stderr,"pepx_processquery: no match for subquery %s\n",subquery);
     return(pepx_reportnomatch(orgquerystring));
     }
   if((cnt=pepx_merge_with_prev_res(finalres, results, acstring, cnt) == 0))
     return(pepx_reportnomatch(orgquerystring));     
   }
//fprintf(stderr,"pepx_processquery last complete subquery OK\n");  
if(i = strlen(query)) // otherwise we're finished
  {
  //fprintf(stderr,"pepx_processquery query not finished, prepare next subquery\n");
  if(strlen(subquery)) // issue last subquery with the longest x-mer
    {
    strncpy(subquery, &querystring[strlen(querystring)-MAXPEPSIZE], MAXPEPSIZE);
    subquery[MAXPEPSIZE] = 0;
    }
  else // query peptide was <= MAXPEPSIZE
   strcpy(subquery,query);
   if(IL_merge)
     cnt = pepx_search(subquery,idxinfoIL);
   else
     cnt = pepx_search(subquery,idxinfo);
  if(cnt==0)
    return(pepx_reportnomatch(orgquerystring));     
  //fprintf(stderr,"last subquery %s gave %d matches\n",subquery,cnt);
  cnt = pepx_merge_with_prev_res(finalres, results, acstring, cnt);
  //fprintf(stderr,"last merge gave %d matches\n",cnt);
  if(cnt == 0)
    return(pepx_reportnomatch(orgquerystring)); 
  }
// last subquery done: finish
// output
pepx_displaymatches(finalres, cnt, orgquerystring, querystring);
return(cnt);  
}

// ---------------- pepx_indexsubseq ---------------------
void pepx_indexsubseq(char* subseq, int varpos)
{
char id[ACLEN];
int i=0, pepsize, sseqcode;
BIN *binptr, *newbinptr, **binarray;

if(varpos) 
  sprintf(id,"%s-%d",currISO,varpos);
else
  strcpy(id,currISO);

pepsize = strlen(subseq);
//if(debug && (pepsize == 6)) fprintf(stderr,"indexing %s (%s)\n",subseq, id);
if(IL_merge)
  // Replace all I/L with J (Was done in pepx_build but needs to be redone for variants)
  for(i=0;i<pepsize;i++)
     if((subseq[i] == 'I') || (subseq[i] == 'L'))
       subseq[i] = 'J';
//fprintf(stderr,"computing sscode for %s\n",subseq);
sseqcode = sscode(subseq);
//fprintf(stderr,"sscode: %X\n",sseqcode);
binarray = toplevel[pepsize];
binptr = binarray[sseqcode];
if(binptr == NULL)
  // 1rst occurence of pep
  {
  binptr = malloc(sizeof(BIN));
  binptr->nextbinptr = NULL;
  binptr->ac_cnt = 1;
  strcpy(binptr->ac[0],id);
  // strcpy(binptr->ac[0],ACpos);  if we later record pos (would double index size)
  binarray[sseqcode] = binptr;
  }
else
  { 
  while(binptr->nextbinptr)
       // follow bins
       binptr = binptr->nextbinptr;
  if(binptr->ac_cnt == BINSIZE)
    {
    // get a new bin
    newbinptr = malloc(sizeof(BIN));
    binptr->nextbinptr = newbinptr;
    binptr = newbinptr;
    binptr->ac_cnt = 0;
    binptr->nextbinptr = NULL;
    //totalbins++;
    }
  i = binptr->ac_cnt++;
  strcpy(binptr->ac[i],id);
  //if(!strcmp(subseq,"AAAA"))  fprintf(stdout,"%s in %s now at %d\n",subseq,currAC,i);
  }
//fprintf(stderr,"pepx_indexsubseq end\n");
}

// ---------------- pepx_indexseq ---------------------
void pepx_indexseq(char* seq, int varcnt)
{
char subseq[MAXPEPSIZE+1], varsubseq[MAXPEPSIZE+1], variant[MAXVARINPEP][MAXPEPSIZE+1]; 
int i, k, varnb, too_many_variants = 0, seqlen, repeat_offset, varindex=1, pepsize;

seqlen = strlen(seq);
if(debug) fprintf(stderr,"indexing %s (%d): %s\n",currISO,seqlen,seq);
for(i=0;i<seqlen-MINPEPSIZE;i++)
   {
   // C-ters
   //printf("indexing at %d/%d\n",i,seqlen);
   //for(pepsize=MINPEPSIZE;pepsize<=MAXPEPSIZE;pepsize++)
   for(pepsize=MINPEPSIZE;pepsize<=MAXPEPSIZE;pepsize++)
      {
      if(seqlen-i < pepsize)
        // subseq too short
        break;
      strncpy(subseq,seq+i,pepsize);
      subseq[pepsize] = 0;
      if(strchr(subseq,'X'))
        // Not indexable
	continue;
      // reference subseq is first in array
      strcpy(variant[0],subseq);
      strcpy(varsubseq,subseq);
      if(debug)
        if(pepsize == 6) fprintf(stderr,"%s at pos %d: %s\n",currAC,i,subseq);	
      if(varcnt)
	// This sequence has variants
	{
        for(k=0;k<pepsize;k++)
          // generate variant subseqs
          {
          strcpy(varsubseq,subseq);
          if(variants[i+k][0])
	    {
	    if(debug)
	      fprintf(stderr,"%s variant at pos %d is %c\n",currISO, i+k,variants[i+k][0]);
	    for(varnb = 0; varnb < MAXVARPERAA; varnb++) 
	       {
	       if(variants[i+k][varnb] == 0)
		        break;
	       if(variants[i+k][varnb] == 'X')
	         // code for missing AA: shift by one and add next
	         {
	         if(k == 0)
	           // replace 1rst by previous
		   *(varsubseq+k) = seq[i-1];
	         else
	           // shift by one and add next
	           {
		   strncpy(varsubseq+k,varsubseq+k+1,pepsize-k);
		   *(varsubseq+pepsize-1) = seq[i+pepsize];
		   }
	         }
	       else
	         *(varsubseq+k) = variants[i+k][varnb];
	       strcpy(variant[varindex++],varsubseq);
 	       //printf("at varindex %d: %s\n",varindex-1, variant[varindex-1]);
	       if(varindex >= MAXVARINPEP)
	         {
	         //fprintf(stderr,"Enough variants for %s in %s\n",subseq,currAC);
	         too_many_variants = 1;
	         varindex_overflow_cnt++;
	         break; // enough
	         }
	       }
	    if(too_many_variants)
	      // Leave pepsize loop
	      {
	      too_many_variants = 0;
	      break;
	      }
	    }
	  }					   
	//printf("varindex loop for %s\n",subseq);
	}
      for(k=0;k<varindex;k++)
        // nonvariant subseq is 0
        {
	strcpy(subseq,variant[k]);
	//printf("varindex loop %d %s\n",k,subseq);
	if(k == 0)
	  // master seq
	  {
	  if(!strstr(seq+i+1,subseq))
            // ensures each ac is indexed only once per sscode
            pepx_indexsubseq(subseq,0); // else  fprintf(stdout,"repeated pep: %s\n",subseq);
	  }
	else if(strlen(subseq) == pepsize && !strstr(masterseq,subseq))
	   // variant peps must not exist elsewhere in master
	   {
	   //if(pepsize == 6) fprintf(stderr,"Variant subseq: %s at pos %d in %s\n",subseq,i+1,currAC); 
     //if(!strcmp(subseq,"TIMTDE")) fprintf(stderr,"masterseq: %s\n",masterseq); 
	   pepx_indexsubseq(subseq,i+1);
	   }
            //else{printf("i:%d k:%d: orgsseq:%s(%d) sseq:%s\n",i,k,variant[0],pepsize,subseq);}
	}
      }
   // Reset varindex for next peptide
   varindex = 1;
   }
}

// ----------------pepx_getFASTAentry ----------------
void pepx_getFASTAentry(char* iso, char* buf, FILE *infile, int* PEFFvarcnt)
{
char *ptr, fastabuf[64];
int seqlen, varcnt, offset;
//fprintf(stderr,"start pepx_getFASTAentry\n");
if(!strncmp(buf,">nxp:",5)) // NetProt peff
  {
  strncpy(iso,buf+8,15);
  *strchr(iso,' ') = 0;
  }
 else // tr fasta
  {
  strncpy(iso,buf+4,15);
  *strchr(iso,'|') = 0;
  strcat(iso,"-1"); // Add an iso id
  }

*buf = 0; // Reset buf to receive sequence
while(fgets(fastabuf,64,infile))
   {
   if(fastabuf[0] == '>') // we reached next header line
      {
      if(!strrchr(buf,'\n')) 
        strcat(buf,"\n"); 
      if(strchr(fastabuf,'\n'))
        offset = -strlen(fastabuf);
      else
        offset = -63;
      //fprintf(stderr,"rewinding with offset %d\n",offset);
      fseek(infile,offset,SEEK_CUR); // rewind line
      return; 
      }
   *strchr(fastabuf,'\n')=0;
   strcat(buf,fastabuf);
   }

 if(!strrchr(buf,'\n'))  
  strcat(buf,"\n");  
 }

// ---------------- pepx_build ---------------------
void pepx_build(char* seqfilename)
{
FILE *in, *varcheck;
char fname[LINELEN], varaa, *ptr, *varptr, varbuf[16], shortenedAC[16];
char buf[MAXHEADSIZE]; // header PEFF line for titin can be over 1.3 Mb
int seqcnt=0, currVarcnt, i, seqlen, maxvar=0, PEFFvarcnt;

if((in=fopen(seqfilename,"r"))==NULL)
  // Sequence file name doesn't exist
 {
 perror(seqfilename);
 exit(2);
 }

if(!ignore_variants)
  {
  fprintf(stdout,"Indexing sequences with NextProt variants, please wait...\n");
  // Check that variant files exist and otherwise build them
  sprintf(fname,"%s/VARIANTS_OK",varfolder);
  if((varcheck=fopen(fname,"r"))==NULL)
    // we must build the files
    {
    if(!pepx_create_variant_files(varfolder))
      exit(55);
    }
  else
    fclose(varcheck);
  }
else
  fprintf(stdout,"Indexing sequences without variants, please wait...\n");
if(IL_merge)
  fprintf(stdout,"Indexing with I/L merged...\n");
pepx_initindexes();
fgets(buf,LINELEN,in); // Skip header line if any
if(strstr(buf,"NX_") || buf[0] == '>')
  {
  // There was no header
  rewind(in);
  if(buf[0] == '>')
    FASTAinput = 1;
  }
else if(buf[0]=='#')
      {
      for(i=0;i<9;i++)
        // PEFF file starts with 10 #-commented description lines
        fgets(buf,MAXHEADSIZE,in);
      FASTAinput = 1;
      }

if(FASTAinput && !ignore_variants)
  {
  fprintf(stderr,"No variants allowed with FASTA input, please relaunch with command flag '--ignore-variants'\n");
  exit(555);
  }

while(fgets(buf,MAXHEADSIZE,in))
    {//fprintf(stderr,"main loop at %d: %s\n",ftell(in),buf);
    if(FASTAinput) // We're parsing the PEFF file
      pepx_getFASTAentry(currISO, buf, in, &PEFFvarcnt);
    else // Parsing from prerelease data
     {
     strncpy(currISO,buf+3,12);
     if(strchr(currISO,'\t'))
      *strchr(currISO,'\t') = 0;
     if(strchr(currISO,' '))
      *strchr(currISO,' ') = 0;
     }

    if(strlen(currISO) > 10 && (!strncmp(currISO,"A0A0",4) || !strncmp(currISO,"A0A1",4)))
      // 10 digit AC, trick it and remember only significant digits in a unique 6-digit length format (there is no real AC starting with P[A-Z]
      // eg: A0A087WTH1 -> PAWTH1
      {
      //fprintf(stderr,"%s -> ",currISO);
      strcpy(shortenedAC,code4tenAA(currISO));
      if(!strcmp(shortenedAC,"XXXXXX"))
      // A new kind of 10 digit AC has appeared
        {
        fprintf(stderr,"Unknown/Untractable AC: %s..skipped\n",currISO);
        continue;
        }
      else  
        strcpy(currISO,shortenedAC);
      //fprintf(stderr,"%s\n",shortenedAC);
      }

    strncpy(currAC,currISO,6);
    currAC[6]=0;
    //if(seqcnt > 30000)      fprintf(stderr,"%s\n",currISO); //debug=TRUE;
    if(FASTAinput)
      strcpy(masterseq,buf);
    else
      strcpy(masterseq,strrchr(buf,'\t')+1);
    while(!isalpha(masterseq[strlen(masterseq)-1])) masterseq[strlen(masterseq)-1] = 0;
    //*strchr(masterseq,'\n')=0; //won't work in PEFF scan for mysterious reason
    seqlen = strlen(masterseq);
    //fprintf(stderr,"%s (%d)\n",currISO,seqlen);
    if(IL_merge)
      // Replace all I/L with J
      for(i=0;i<seqlen;i++)
         if((masterseq[i] == 'I') || (masterseq[i] == 'L'))
            masterseq[i] = 'J';

    if(!strcmp(currAC,"Pxxxxx"))
    //if(!strcmp(currAC,"PBJ2A2"))
      {
      debug=TRUE;
      fprintf(stderr,"input buffer: %s...",buf);
      }
    else			
      debug=FALSE;		  		      
    //debug=TRUE;		  		      
    // get variants for current entry
    if(!ignore_variants)
      {
       if(!FASTAinput) 
         currVarcnt = pepx_build_varindex(currISO);
       else
        currVarcnt = PEFFvarcnt; // varindex already built while processing PEFF entry in pepx_getPEFFvariants
      }
    else
       currVarcnt = 0;
    totalvarcnt += currVarcnt;
    //debug=TRUE;
    pepx_indexseq(masterseq,currVarcnt);
    seqcnt++;
    if((seqcnt % 1000) == 0)
      fprintf(stderr,"Processing seq %d...\n",seqcnt);
    }
fprintf(stderr,"%d sequences indexed (%d simple variants, %d/%d/%d 1-miss variants)\n",seqcnt,totalvarcnt,varmisscnt1, varmisscnt2, varmisscnt3);
fprintf(stderr,"%d varindex overflows\n",varindex_overflow_cnt);
if(!debug)
  pepx_saveall();
}

// ---------------- printHelp ---------------------
void  printHelp(char *mode)
{
if(!strcmp(mode,"ARGS"))
  {
  fprintf(stderr,"\nPepx possible arguments are:\n\n");
  fprintf(stderr,"--search (short=-s) to perform a peptide search\n");
  fprintf(stderr,"--build (short=-b) to build indexes (requires the isoform-file name as mandatory argument)\n");
  fprintf(stderr,"--version (short=-v) to show current version\n");
  fprintf(stderr,"--help (short=-h) to show this help\n");
  fprintf(stderr,"--index-folder (short=-x) to specify an index folder (default is .)\n");
  fprintf(stderr,"--variant-folder (short=-w) to specify a folder for json variants (required for build command when ignore-variants flag is not set)\n");
  fprintf(stderr,"--json output in json format (incompatible with -noiso mode)\n");
  fprintf(stderr,"--peptide-file (short=-p) a file with peptides to search (1 peptide/line, if not provided peptides will be read from stdin)\n");
  fprintf(stderr,"--ignore-variants to build indexes not considering variants\n");
  fprintf(stderr,"--IL to build indexes merging I and L\n");
  fprintf(stderr,"--noiso (short=-n) to output search results at the entry level\n");
  }
fprintf(stderr,"\nCurrent limitation:\n\n");
fprintf(stderr,"- poly-AA stretches > %d cannot be found\n",MAXPEPSIZE);
fprintf(stderr,"- only snp-style (1 AA for 1 other AA), and 1-AA-miss variants are accounted \n");
fprintf(stderr,"- only 32 variant accounted within a given x-mer\n");
fprintf(stderr,"- only 1 joker (X) allowed in a given x-mer\n"); 
}

// ---------------- main ---------------------
int main(int argc, char **argv)
{
char *ptr, *peptoken, seqfname[LINELEN]="", command[LINELEN], querystring[8192];
char query[8192], pepfname[LINELEN]="";
int i, c;
int option_index = 0; // getopt_long stores the option index here.
FILE* inputstream = stdin;

static struct option long_options[] = {
               {"ignore-variants",   no_argument,       &ignore_variants, 1},
               {"IL",   no_argument,       &IL_merge, 1},
               {"search",  no_argument,       0, 's'},
               {"noiso",   no_argument,       0, 'n'},
               {"help",   no_argument,       0, 'h'},
               {"version",   no_argument,       0, 'v'},
               {"json",   no_argument,       &json, 1},
               {"build", required_argument, 0, 'b'},
               {"peptide-file", required_argument, 0, 'p'},
               {"variant-folder",  required_argument, 0, 'w'},
               {"index-folder",required_argument, 0, 'x'},
               {0, 0, 0, 0}
             };
// Further improvements: Rewrite the code in scala to take advantage of powerful hashmap-like data structures
// Consider indexing non-snp variants (alleles...)

if(argc < 2)
 {
 if((ptr=getenv("QUERY_STRING")) == NULL)
   {
   fprintf(stderr,"Usage: pepx [--build or --search] [filename or peptide or INTERACTIVE (if search)]\n");
   fprintf(stderr,"if no peptide is given as second arg, then they are read from stdin\n");
   fprintf(stderr,"type pepx -h for full documentation\n\n");
   exit(1);
   }
 // for web cgi
 strcpy(envstring,ptr);
 fprintf(stdout,"Content-type: text/html\n\n");
 if(ptr=strstr(envstring,"pep="))
   {
   strcpy(command,"search");
   strcpy(querystring,ptr+4);
   strcpy(linesep,"<br>");
   if(ptr=strstr(envstring,"=noiso"))
      // output matches at the entry level
      strcpy(matchmode,"ACONLY");
   //if(ptr=strstr(envstring,"=IL"))
   IL_merge = 1; // IL-merged indexes is the default mode for the web
   if(ptr=strstr(envstring,"=json"))
      // output will be in json
      json = 1;
   
   }
 else
   {
   fprintf(stdout,"envstring: %s\n",envstring);
   fprintf(stdout,"No valid arguments...\n");
   exit(0);
   }
 }
else
  // non-web usage: parse command arguments
  {
  //while((c = getopt_long (argc, argv, "snhvb:p:f:w:r:x:", long_options, &option_index)) != -1)
  while((c = getopt_long (argc, argv, "snhvb:p:f:w:x:", long_options, &option_index)) != -1)
        {
       switch (c)
             {
             case 0:
               break;
     
             case 'n':
               strcpy(matchmode,"ACONLY"); // For search command only
               break;
     
             case 'v':
	       printf ("Current version is `%s'\n", version);
               exit(0);
               
             case 'h':
               printHelp("ARGS");
               exit(0);
     
             case 's':
	       strcpy(command,"search");
               break;
               
             case 'b':
               strcpy(command,"build");
	       strcpy(seqfname,optarg);
               break;
     
             case 'p':
	       strcpy(pepfname,optarg);
	       if((inputstream=fopen(pepfname,"r"))==NULL)
	         {
		 perror(pepfname);
	         exit(6);
		 }
               break;
     
             case 'w':
               strcpy(varfolder,optarg);
               break;
     
             case 'x':
               strcpy(indexfolder, optarg);
               break;
     
             default:
	       printf("%c: Unknown option...\n", c);
               abort();
             }

       }
  }

if(!strcmp(command,"build"))
  {
  if(strlen(seqfname) == 0)
    {
    fprintf(stderr,"'build' requires a sequence file for argument\n");
    exit(2);
    }
  if ((ignore_variants == 0) && strlen(varfolder) == 0)  
    {
    fprintf(stderr,"'build' with variants requires a variant folder for argument (option -w)\n");
    exit(2);
    }
  start = clock();
  //strcpy(seqfname,argv[2]);
  printf("Building with seqfile %s, indexes at %s, variants at %s\n",seqfname,indexfolder,ignore_variants?"ignored":varfolder);
  //exit(0);
  pepx_build(seqfname);
  fprintf(stderr, "Duration: %d\n", clock() - start);
  }
else if(!strncmp(command,"search",6))
  {
  if(optind < argc)
    strcpy(querystring,argv[optind]);
  else if(strlen(envstring) == 0)
    strcpy(querystring,"BATCH");
  if(!json)
     printf("Searching for %s, indexes at %s\n",querystring,indexfolder);
  //exit(0);
  pepx_loadall();
  if(!strcmp(querystring,"INTERACTIVE") || !strcmp(querystring,"BATCH"))
    {
    strcpy(outputmode,querystring);
    printHelp("");
    do
     {
     if(!strcmp(querystring,"INTERACTIVE"))
       printf("\nEnter or paste peptide query (spaces and digits will be skipped)? ");
     // get user input
     if(!fgets(query,8192,inputstream))
       break;
     if(ptr=strrchr(query,'\n'))
       *ptr=0;
     if(strlen(query) < 2)
       exit(0);
     pepx_processquery(query);
     }
    while(TRUE);
    }
  else
     // Just process given query(ies)
     {
     if(json)
       {
       fprintf(stdout,"{\n"); // opens json stream
       pepx_json_header("params",querystring,IL_merge);
       fprintf(stdout,"\"peptideMatches\":[");
       }

     peptoken = strtok(querystring, ",");  // there may be several comma-separated peptides 
     while( peptoken != NULL )
          {
          pepx_processquery(peptoken);
          peptoken = strtok(NULL, ",");
          if(json && (peptoken != NULL)) // prepare next entryematches
	    fprintf(stdout,",\n");
	  }
     if(json) fprintf(stdout,"]\n}\n"); // close  peptidematches  and json object	
     }
  }
else
  {
  fprintf(stderr,"command arg must be either 'build' or 'search'\n");
  exit(4);
  }
return(0);
}


