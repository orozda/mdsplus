#include        <stdarg.h>
#include        <stdio.h>
#include        <stdlib.h>
#if defined(vms)
#include        <lib$routines.h>
#include        <time.h>
#include        <unistd.h>
#elif defined(_WIN32)
#include        <io.h>
#define isatty(a) _isatty(a)
#elif defined(vxWorks)
#include 	<time.h>
#else
#include        <sys/time.h>
#endif

#define CREATE_STS_TEXT
#include        "mdsdcldef.h"
#include        "clidef.h"
#include        "tcldef.h"
#include        "ccldef.h"
#include        "librtl_messages.h"
#include        "tdimessages.h"
#include        "treeshr.h"
#include        "servershr.h"

/**********************************************************************
* MDSMSG.C --
*
* Routine for printing status information ...
*
* History:
*  06-Apr-2001  TRG  Start/end each MdsMsg line with '\r', esp for IDL.
*  05-Jan-1998  TRG  Create.
*
************************************************************************/


		/*========================================================
		 * "Define"s and structure definitions ...
		 *=======================================================*/
#define ALREADY_DISPLAYED  0x80000000


		/*========================================================
		 * Function prototypes ...
		 *=======================================================*/
extern void StrCopyDx();
extern int MDSprintf( char *fmt , ... );
extern int MDSfprintf( FILE *fp , char *fmt , ... );


		/*========================================================
		 * Static variables ...
		 *=======================================================*/
static int   mdsmsgFlag = 1;		/* 1 for longer "status" string	*/



	/*****************************************************************
	 * getFacility:
	 *****************************************************************/
static int   getFacility(	/* Return: num entries in stsText[]	*/
    int   sts			/* <r> sts value			*/
   ,char  **facilityText	/* <w> name of facility: c-string	*/
   ,struct stsText  **stsText	/* <w> addr of facility's "stsText[]"	*/
   )
   {
    int   facility;
    int   max;

    facility = (sts & ~ALREADY_DISPLAYED) >> 16;
    if (facility == MDSDCL_FACILITY)
       {
        *stsText = mdsdcl_stsText;	/* point to array		*/
        *facilityText = "MDSDCL_FACILITY";
        max = sizeof(mdsdcl_stsText)/sizeof(mdsdcl_stsText[0]);
       }
    else if (facility == CLI_FACILITY)
       {
        *stsText = cli_stsText;		/* point to array		*/
        *facilityText = "CLI_FACILITY";
        max = sizeof(cli_stsText)/sizeof(cli_stsText[0]);
       }
    else if (facility == CCL_FACILITY)
       {
        *stsText = ccl_stsText;		/* point to array		*/
        *facilityText = "CCL_FACILITY";
        max = sizeof(ccl_stsText)/sizeof(ccl_stsText[0]);
       }
    else if (facility == SERVERSHR_FACILITY)
       {
        *stsText = servershr_stsText;		/* point to array		*/
        *facilityText = "SERVERSHR_FACILITY";
        max = sizeof(servershr_stsText)/sizeof(servershr_stsText[0]);
       }
    else if (facility == TDI_FACILITY)
       {
        *stsText = tdi_stsText;		/* point to array		*/
        *facilityText = "TDI_FACILITY";
        max = sizeof(tdi_stsText)/sizeof(tdi_stsText[0]);
       }
    else if (facility == TCL_FACILITY)
       {
        *stsText = tcl_stsText;		/* point to array		*/
        *facilityText = "TCL_FACILITY";
        max = sizeof(tcl_stsText)/sizeof(tcl_stsText[0]);
       }
    else if (facility == TREESHR_FACILITY)
       {
        *stsText = treeshr_stsText;	/* point to array		*/
        *facilityText = "TREESHR_FACILITY";
        max = sizeof(treeshr_stsText)/sizeof(treeshr_stsText[0]);
       }
    else if (facility == LIB_FACILITY)
       {
        *stsText = librtl_stsText;	/* point to array		*/
        *facilityText = "LIB_FACILITY";
        max = sizeof(librtl_stsText)/sizeof(librtl_stsText[0]);
       }
    else if (facility == STRMDS_FACILITY)
       {
        *stsText = strMds_stsText;	/* point to array		*/
        *facilityText = "STRMDS_FACILITY";
        max = sizeof(strMds_stsText)/sizeof(strMds_stsText[0]);
       }
    else if (facility == SS_FACILITY)
       {
        *stsText = ss_stsText;	/* point to array		*/
        *facilityText = "SS_FACILITY";
        max = sizeof(ss_stsText)/sizeof(ss_stsText[0]);
       }
    else
       {
        *stsText = 0;
        *facilityText = "???";
        max = 0;
       }

    return(max);
   }



	/*****************************************************************
	 * MdsGetMsg:
	 *****************************************************************/
char  *MdsGetMsg(	/* Return: addr of "status" string	*/
    int   sts			/* <r> sts value			*/
   )
   {
    int   i;
    int   max;
    char  *facilityText;
    struct stsText  *stsText;
    static char  text[72];
    static DESCRIPTOR(dsc_text,text);

    max = getFacility(sts,&facilityText,&stsText);

    for (i=0 ; i<max ; i++)
       {
        if (sts == stsText[i].stsL_num)
          {
           str_copy_dx(&dsc_text,stsText[i].stsA_name);
           if (mdsmsgFlag & 1)
              {		/* if longer messages ...		*/
               str_append(&dsc_text," : ");
               str_append(&dsc_text,stsText[i].stsA_text);
              }
           return(text);
          }
       }

    sprintf(text,"%s %s : %08X",facilityText,(sts&1)?"Success":"Error",sts);
    return(text);
   }


void MdsGetMsgDsc(int status, struct descriptor *out)
{
	static struct descriptor tmp = {0,DTYPE_T,CLASS_S,0};
	tmp.dscA_pointer = MdsGetMsg(status);
	tmp.dscW_length = strlen(tmp.dscA_pointer);
	StrCopyDx(out,&tmp);
	return;
}

	/****************************************************************
	 * setMdsMsgFlag:
	 ****************************************************************/
int   setMdsMsgFlag(		/* Return: old setting of flag		*/
    int   new			/* <r> new setting for flag		*/
   )
   {
    int   old;

    old = mdsmsgFlag;
    mdsmsgFlag = new;
    return(old);
   }



	/*****************************************************************
	 * MdsMsg:
	 *****************************************************************/
int   MdsMsg(			/* Return: sts provided by user		*/
    int   sts			/* <r> status code			*/
   ,char  fmt[]			/* <r> format statement			*/
   , ...			/* <r:opt> arguments to fmt[]		*/
   )
   {
    int   k;
    int   write2stdout;
    char  text[256];
    va_list  ap;		/* arg ptr				*/

    write2stdout = isatty(fileno(stdout)) ^ isatty(fileno(stderr));
    if ((sts & ALREADY_DISPLAYED) && (sts != -1))
        return(sts);

    sprintf(text,"\r%s: ",pgmname());
    k = strlen(text);
    if (fmt)
       {
        va_start(ap,fmt);	/* initialize "ap"			*/
        vsprintf(text+k,fmt,ap);
        if (sts)
           {
            MDSfprintf(stderr,"%s\n\r    sts=%s\n\n\r",text,MdsGetMsg(sts));
            if (write2stdout)
                MDSfprintf(stdout,"%s\n\r    sts=%s\n\n\r",
                    text,MdsGetMsg(sts));
           }
        else
           {
            MDSfprintf(stderr,"%s\n\r",text);
            if (write2stdout)
                MDSfprintf(stdout,"%s\n\r",text);
           }
       }
    else
       {
        MDSfprintf(stderr,"%s:  sts=%s\n\r",text,MdsGetMsg(sts));
        if (write2stdout)
            MDSfprintf(stdout,"%s:  sts=%s\n\r",text,MdsGetMsg(sts));
       }

    return(sts | ALREADY_DISPLAYED);
   }

#ifdef MAIN
void  main()
   {
    MdsMsg(MDSDCL_STS_SUCCESS,0);  MDSprintf("\n");
    MdsMsg(CLI_STS_PRESENT,0);  MDSprintf("\n");
    MdsMsg(CCL_STS_SUCCESS,0);  MDSprintf("\n");
    MdsMsg(TCL_STS_SUCCESS,0);  MDSprintf("\n");
    MdsMsg(TreeNORMAL,0);  MDSprintf("\n");
    MdsMsg(LibNOTFOU,0);  MDSprintf("\n");
    MdsMsg(StrMATCH,0);  MDSprintf("\n");
    MdsMsg(SsINTOVF,0);  MDSprintf("\n");
   }
#endif
