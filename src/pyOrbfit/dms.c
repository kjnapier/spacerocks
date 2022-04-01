/* 	$Id: dms.c,v 2.0 2001/09/14 19:05:06 garyb Exp $	 */

#ifndef lint
static char vcid[] = "$Id: dms.c,v 2.0 2001/09/14 19:05:06 garyb Exp $";
#endif /* lint */
/* convert a string to decimal degrees */
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>

double 
dmsdeg(char *string)
{
	float args[3];
	char substring[10];
	int i,iarg,j,sgn;
	for (iarg=0;iarg<3;args[iarg++]=0.);

	/*strip off leading whitespace and look for - sign*/
	for (i=0; isspace( (int) string[i]); i++);
	if (string[i]=='-') {
		sgn = -1;
		i++;
	}
	else {
		sgn = 1;
	}

	for (iarg=0; iarg<3; iarg++) {
	  substring[0] = 0;
	  j = sscanf(string+i, "%[ :]",substring);
	  if (j==EOF || *(string+i)==0) break;  /*failing to get EOF at end*/
	  i += strlen(substring);
	  substring[0] = 0;
	  j = sscanf(string+i,"%[-+.0-9]",substring);
	  if (j==EOF || *(string+i)==0) break;
	  else if (j!=1) {
	    fprintf(stderr,"Error interpeting dms ->%s<-\n",string);
	    return(0.);
	  }
	  args[iarg] = atof(substring);
	  i += strlen(substring);
	}
	args[1] += args[2]/60.;
	if (args[0]<0) args[1] = - args[1];
	args[0] += args[1]/60.;
	return( args[0]*sgn );
}

/* slight variation for hours:*/
double 
hmsdeg(char *string)
{
	return(15. * dmsdeg(string));
}

void
degdms(double degr,
       char *outbuff)
{
  int  ideg,imin,sign=1;
  double sec;

  if (degr<0.) {
    sign=-1;
    degr *= -1.;
  }

  ideg = floor(degr);
  degr = (degr-ideg)*60.;
  imin = floor(degr);
  sec = (degr-imin)*60.;

  if (sign>0) 
    if (ideg<100)
      sprintf(outbuff," +%02d:%02d:%06.3f",ideg,imin,sec);
    else
      sprintf(outbuff,"+%03d:%02d:%06.3f",ideg,imin,sec);
  else
    if (ideg<100)
      sprintf(outbuff," -%02d:%02d:%06.3f",ideg,imin,sec);
    else
      sprintf(outbuff,"-%03d:%02d:%06.3f",ideg,imin,sec);
  return;
}
void
deghms(double degr,
       char *outbuff)
{
  int  ideg,imin,sign=1;
  double sec;

  if (degr<0) {
    sign=-1;
    degr *= -1.;
  }

  degr /= 15.;
  ideg = floor(degr);
  degr = (degr-ideg)*60.;
  imin = floor(degr);
  sec = (degr-imin)*60.;

  if (sign>0) 
    sprintf(outbuff," %02d:%02d:%07.4f",ideg,imin,sec);
  else
    sprintf(outbuff,"-%02d:%02d:%07.4f",ideg,imin,sec);

  return;
}

