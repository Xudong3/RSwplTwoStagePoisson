


double FouOrdPiInternal(int pos1, int pos2, int pos3, int pos4, int sc1, int sc2,int sc3, int sc4, double fss[], double n2infor[],double N2){

   if (sc1==sc2 & sc2==sc3 & sc3==sc4){
      if (pos1==pos2 & pos2==pos3 & pos3==pos4){
	return n2infor[sc1]/N2*fss[sc1];
      } else if ((pos1==pos2) & !(pos2==pos3) & pos3==pos4){
	 return n2infor[sc1]/N2* (n2infor[sc1]-1)/(N2-1)*fss[sc1];
      }else if ((pos1==pos3) & !(pos3==pos2) & pos2==pos4){
	 return n2infor[sc1]/N2*(n2infor[sc1]-1)/(N2-1)*fss[sc1];
      } else if ((pos1==pos4) & !(pos4==pos2) & pos2==pos3){
	 return n2infor[sc1]/N2*(n2infor[sc1]-1)/(N2-1)*fss[sc1];
      } else if (pos1==pos2 & pos2==pos3 & !(pos3==pos4)){
	 return n2infor[sc1]/N2* (n2infor[sc1]-1)/(N2-1)*fss[sc1];
      } else if ((pos1==pos2) & (pos2==pos4) & !(pos4==pos3)){
	 return n2infor[sc1]/N2*(n2infor[sc1]-1)/(N2-1)*fss[sc1];
      }else if ((pos1==pos3) & (pos3==pos4) & !(pos4==pos2)){
	 return n2infor[sc1]/N2*(n2infor[sc1]-1)/(N2-1)*fss[sc1];
      } else if (pos2==pos3 & pos3==pos4 & !(pos4==pos1)){
	 return n2infor[sc1]/N2*(n2infor[sc1]-1)/(N2-1)*fss[sc1];
      } else if ((pos1==pos2) & !(pos2==pos3) & !(pos4==pos3)  &!(pos4==pos1)){
	 return n2infor[sc1]/N2*(n2infor[sc1]-1)/(N2-1)*(n2infor[sc1]-2)/(N2-2) *fss[sc1];
      } else if ((pos1==pos3) & !(pos2==pos1) & !(pos4==pos2)  &!(pos4==pos1)){
	 return n2infor[sc1]/N2*(n2infor[sc1]-1)/(N2-1)*(n2infor[sc1]-2)/(N2-2)*fss[sc1];
      } else if ((pos1==pos4) & !(pos2==pos1) & !(pos3==pos2)  &!(pos3==pos1)){
	return  n2infor[sc1]/N2*(n2infor[sc1]-1)/(N2-1)*(n2infor[sc1]-2)/(N2-2)*fss[sc1];
      } else if ((pos2==pos3) & !(pos1==pos2) & !(pos4==pos1)  &!(pos4==pos2)){
	 return n2infor[sc1]/N2*(n2infor[sc1]-1)/(N2-1)*(n2infor[sc1]-2)/(N2-2)*fss[sc1];
      }else if ( (pos2==pos4) & !(pos1==pos2) & !(pos3==pos1)  &!(pos3==pos2)){
         return  n2infor[sc1]/N2*(n2infor[sc1]-1)/(N2-1)*(n2infor[sc1]-2)/(N2-2)*fss[sc1];
      } else if ((pos3==pos4) & !(pos1==pos3) & !(pos2==pos1)  &!(pos2==pos3)){
	 return n2infor[sc1]/N2*(n2infor[sc1]-1)/(N2-1)*(n2infor[sc1]-2)/(N2-2)*fss[sc1];
      } else {
	 return n2infor[sc1]/N2*(n2infor[sc1]-1)/(N2-1)*(n2infor[sc1]-2)/(N2-2)*(n2infor[sc1]-3)/(N2-3)*fss[sc1];
      }
   }else if ((sc1==sc2) & !(sc2==sc3) & sc3==sc4){
      if (!(pos1==pos2) & !(pos3==pos4) ){
         return n2infor[sc1]/N2*(n2infor[sc1]-1)/(N2-1)*n2infor[sc3]/N2*(n2infor[sc3]-1)/(N2-1)*fss[sc1]*fss[sc3];
      }else if ((pos1==pos2) & !(pos3==pos4)) {
	 return n2infor[sc1]/N2*n2infor[sc3]/N2* (n2infor[sc3]-1)/(N2-1)*fss[sc1]*fss[sc3];
      }else if (!(pos1==pos2) & (pos3==pos4)){
	 return n2infor[sc1]/N2* (n2infor[sc1]-1)/(N2-1)* n2infor[sc3]/N2*fss[sc1]*fss[sc3];
      }else if (pos1==pos2 & pos3==pos4){
	 return n2infor[sc1]/N2* n2infor[sc3]/N2*fss[sc1]*fss[sc3];
      }
   }else if ( (sc1==sc3) & !(sc3==sc2) & sc2==sc4){
      if (!(pos1==pos3) & !(pos2==pos4)){
	 return n2infor[sc1]/N2*(n2infor[sc1]-1)/(N2-1)* n2infor[sc2]/N2*(n2infor[sc2]-1)/(N2-1)*fss[sc1]*fss[sc2];
      }else if ((pos1==pos3) & !(pos2==pos4)){
	 return n2infor[sc1]/N2* n2infor[sc2]/N2*(n2infor[sc2]-1)/(N2-1) *fss[sc1]*fss[sc2];
      } else if ( !(pos1==pos3) & (pos2==pos4)){
	return  n2infor[sc1]/N2*(n2infor[sc1]-1)/(N2-1)* n2infor[sc2]/N2*fss[sc1]*fss[sc2];
      }  else if (pos1==pos3 & pos2==pos4){
	return  n2infor[sc1]/N2* n2infor[sc2]/N2*fss[sc1]*fss[sc2];
      }
   }else if ((sc1==sc4) & !(sc4==sc2) & sc2==sc3){
      if (!(pos1==pos4) & !(pos2==pos3)){
	 return n2infor[sc1]/N2*(n2infor[sc1]-1)/(N2-1)*n2infor[sc2]/N2*(n2infor[sc2]-1)/(N2-1)*fss[sc1]*fss[sc2];
      } else if ((pos1==pos4) & !(pos2==pos3)){
	 return n2infor[sc1]/N2*n2infor[sc2]/N2*(n2infor[sc2]-1)/(N2-1)*fss[sc1]*fss[sc2];
      } else if (!(pos1==pos4) & (pos2==pos3) ){
	 return n2infor[sc1]/N2*(n2infor[sc1]-1)/(N2-1)* n2infor[sc2]/N2*fss[sc1]*fss[sc2];
      } else if (pos1==pos4 & pos2==pos3){
	 return n2infor[sc1]/N2* n2infor[sc2]/N2*fss[sc1]*fss[sc2];
      }
   }else if (sc1==sc2 & sc2==sc3 & !(sc3==sc4)){
      if (!(pos1==pos2)& !(pos2==pos3) & !(pos3==pos1)){
	 return n2infor[sc1]/N2*(n2infor[sc1]-1)/(N2-1)*(n2infor[sc1]-2)/(N2-2)* n2infor[sc4]/N2*fss[sc1]*fss[sc4];
      }  else if ( (pos1==pos2) & !(pos2==pos3)){
	return  n2infor[sc1]/N2*(n2infor[sc1]-1)/(N2-1)*n2infor[sc4]/N2*fss[sc1]*fss[sc4];
      }  else if ( (pos1==pos3) & !(pos3==pos2)){
	return  n2infor[sc1]/N2*(n2infor[sc1]-1)/(N2-1)* n2infor[sc4]/N2*fss[sc1]*fss[sc4];
      }  else if ( (pos2==pos3) & !(pos3==pos1)){
	return  n2infor[sc1]/N2*(n2infor[sc1]-1)/(N2-1)*n2infor[sc4]/N2*fss[sc1]*fss[sc4];
      }  else if ( pos2==pos3 & pos3==pos1){
	return  n2infor[sc1]/N2*n2infor[sc4]/N2*fss[sc1]*fss[sc4];
      }
   }else if (sc1==sc2 & sc2==sc4 & !(sc4==sc3)){
      if (!(pos1==pos2) & !(pos2==pos4) & !(pos4==pos1) ){
	 return n2infor[sc1]/N2*(n2infor[sc1]-1)/(N2-1)*(n2infor[sc1]-2)/(N2-2)* n2infor[sc3]/N2*fss[sc1]*fss[sc3];
      } else if ((pos1==pos2) & !(pos2==pos4)){
	return n2infor[sc1]/N2*(n2infor[sc1]-1)/(N2-1)*n2infor[sc3]/N2*fss[sc1]*fss[sc3];  /* ?? */
      } else if ((pos1==pos4) & !(pos4==pos2)){
        return  n2infor[sc1]/N2*(n2infor[sc1]-1)/(N2-1)*
	   n2infor[sc3]/N2*fss[sc1]*fss[sc3];
      } else if ((pos2==pos4) & !(pos4==pos1)){
        return  n2infor[sc1]/N2*(n2infor[sc1]-1)/(N2-1)*
	   n2infor[sc3]/N2*fss[sc1]*fss[sc3];
      } else if (pos1==pos2 & pos2==pos4){
        return  n2infor[sc1]/N2*
	   n2infor[sc3]/N2*fss[sc1]*fss[sc3];
      }
   }else if (sc1==sc3 & sc3==sc4 & !(sc4==sc2)){
      if (!(pos1==pos3) & !(pos3==pos4) & !(pos4==pos1)){
         return n2infor[sc1]/N2*(n2infor[sc1]-1)/(N2-1)*(n2infor[sc1]-2)/(N2-2)*
	   n2infor[sc2]/N2*fss[sc1]*fss[sc2];
      } else if ((pos1==pos3) & !(pos3==pos4)){
         return n2infor[sc1]/N2*(n2infor[sc1]-1)/(N2-1)*
	   n2infor[sc2]/N2*fss[sc1]*fss[sc2];
      } else if ((pos1==pos4) & !(pos4==pos3)){
         return n2infor[sc1]/N2*(n2infor[sc1]-1)/(N2-1)*
	   n2infor[sc2]/N2*fss[sc1]*fss[sc2];
      } else if ((pos3==pos4) & !(pos4==pos1)){
        return  n2infor[sc1]/N2*(n2infor[sc1]-1)/(N2-1)*
	   n2infor[sc2]/N2*fss[sc1]*fss[sc2];
      }else if (pos3==pos4 & pos4==pos1 ){
        return  n2infor[sc1]/N2*
	   n2infor[sc2]/N2*fss[sc1]*fss[sc2];
      }
   }else if (sc2==sc3 & sc3==sc4 & !(sc4==sc1)){
      if (!(pos2==pos3) & !(pos3==pos4) & !(pos4==pos2)){
        return  n2infor[sc2]/N2*(n2infor[sc2]-1)/(N2-1)*(n2infor[sc2]-2)/(N2-2)*
	   n2infor[sc1]/N2*fss[sc2]*fss[sc1];
      }else if ((pos2==pos3) & !(pos3==pos4)){
        return  n2infor[sc2]/N2*(n2infor[sc2]-1)/(N2-1)*
	   n2infor[sc1]/N2*fss[sc2]*fss[sc1];
      } else if ((pos2==pos4) & !(pos4==pos3)){
        return  n2infor[sc2]/N2*(n2infor[sc2]-1)/(N2-1)*
	   n2infor[sc1]/N2*fss[sc2]*fss[sc1];
      } else if ((pos3==pos4) & !(pos4==pos2)){
        return  n2infor[sc2]/N2*(n2infor[sc2]-1)/(N2-1)*
	   n2infor[sc1]/N2*fss[sc2]*fss[sc1];
      }else if (pos3==pos4 & pos4==pos2){
         return n2infor[sc2]/N2*
	   n2infor[sc1]/N2*fss[sc2]*fss[sc1];
      }
   }else if ((sc1==sc2) & !(sc2==sc3) & !(sc4==sc3)  &!(sc4==sc1)){
      if (!(pos1==pos2)){
        return  n2infor[sc1]/N2*(n2infor[sc1]-1)/(N2-1)*
	   n2infor[sc3]/N2*n2infor[sc4]/N2*fss[sc1]*fss[sc3]*fss[sc4];
      } else if  (pos1==pos2){
         return n2infor[sc1]/N2*
	   n2infor[sc3]/N2*n2infor[sc4]/N2*fss[sc1]*fss[sc3]*fss[sc4];
      }
   }else if ((sc1==sc3) & !(sc2==sc1) & !(sc4==sc2)  &!(sc4==sc1)){
      if (!(pos1==pos3)){
        return  n2infor[sc1]/N2*(n2infor[sc1]-1)/(N2-1)*
	   n2infor[sc2]/N2*n2infor[sc4]/N2*fss[sc1]*fss[sc2]*fss[sc4];
      } else if (pos1==pos3){
        return n2infor[sc1]/N2*
	  n2infor[sc2]/N2 *n2infor[sc4]/N2*fss[sc1]*fss[sc2]*fss[sc4];
      }
   }else if ((sc1==sc4) & !(sc2==sc1) & !(sc3==sc2)  &!(sc3==sc1)){
      if (!(pos1==pos4)){
        return  n2infor[sc1]/N2*(n2infor[sc1]-1)/(N2-1)*
	   n2infor[sc2]/N2*n2infor[sc3]/N2*fss[sc1]*fss[sc2]*fss[sc3];
      }else if (pos1==pos4){
         return n2infor[sc1]/N2*
	   n2infor[sc2]/N2*n2infor[sc3]/N2*fss[sc1]*fss[sc2]*fss[sc3];
      }
   }else if ((sc2==sc3) & !(sc1==sc2) & !(sc4==sc1)  &!(sc4==sc2)){
      if (!(pos2==pos3)){
         return n2infor[sc2]/N2*(n2infor[sc2]-1)/(N2-1)*
	   n2infor[sc1]/N2*n2infor[sc4]/N2*fss[sc2]*fss[sc1]*fss[sc4];
      }else if (pos2==pos3){
         return n2infor[sc2]/N2*
	   n2infor[sc1]/N2*n2infor[sc4]/N2*fss[sc2]*fss[sc1]*fss[sc4];
      }
   }else if ((sc2==sc4) & !(sc1==sc2) & !(sc3==sc1)  &!(sc3==sc2)){
      if (!(pos2==pos4)){
         return n2infor[sc2]/N2*(n2infor[sc2]-1)/(N2-1)*
	   n2infor[sc1]/N2*n2infor[sc3]/N2*fss[sc2]*fss[sc1]*fss[sc3];
      }else if (pos2==pos4){
         return n2infor[sc2]/N2*
	   n2infor[sc1]/N2*n2infor[sc3]/N2*fss[sc2]*fss[sc1]*fss[sc3];
      }
   }else if ((sc3==sc4) & !(sc1==sc3) & !(sc2==sc1)  &!(sc2==sc3)){
      if (!(pos3==pos4)){
         return n2infor[sc3]/N2*(n2infor[sc3]-1)/(N2-1)*
	   n2infor[sc1]/N2*n2infor[sc2]/N2*fss[sc3]*fss[sc1]*fss[sc2];
      }else if (pos3==pos4){
         return n2infor[sc3]/N2*
	   n2infor[sc1]/N2*n2infor[sc2]/N2*fss[sc3]*fss[sc1]*fss[sc2];
   }
      }else{
      return n2infor[sc1]/N2*n2infor[sc2]/N2*
	n2infor[sc3]/N2*n2infor[sc4]/N2*fss[sc1]*fss[sc2]*fss[sc3]*fss[sc4];
   }
   return -1.0;
}

double SecOrdPiInternal(int pos1, int pos2, int sc1, int sc2, double fss[], double n2infor[], double N2){
	
  if (pos1==pos2 & sc1==sc2){
      return n2infor[sc1]/N2*fss[sc1];
   }else if (!(pos1==pos2) & (sc1==sc2)){
      return n2infor[sc1]/N2*(n2infor[sc1]-1)/(N2-1)*fss[sc1];
   }else if ((pos1==pos2) & !(sc1==sc2)){
      return 0;  /* can't happen */
   }else if (!(pos1==pos2) & !(sc1==sc2)){
      return n2infor[sc1]/N2*n2infor[sc2]/N2*fss[sc1]*fss[sc2];
   }
   return -1; /* can't happen */
}

void SecOrdPi(int pos1[], int pos2[], int sc1[], int sc2[], double  fss[], double n2infor[], double *N2, int *len, double answer[]){
	int i;
		 
    for(i=0;i<*len; i++){
      answer[i] = SecOrdPiInternal(pos1[i], pos2[i], 
				 sc1[i]-1, sc2[i]-1, fss,
				 n2infor, *N2);
     }
    return;
}	

void FourOrdPi( int pos1[], int pos2[], int pos3[], int pos4[], 
                 int sc1[], int sc2[], int sc3[], int sc4[], 
                 double fss[], double n2infor[], double *N2,
				 int *len, double answer[]){
  int i;

  for(i=0;i<*len; i++){
    answer[i] = FouOrdPiInternal(pos1[i], pos2[i], pos3[i], pos4[i], 
				 sc1[i]-1, sc2[i]-1, sc3[i]-1,sc4[i]-1,
				fss,  n2infor, *N2);

  }
  return;
}
