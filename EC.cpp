
#include <iostream>
#include <vector>
#include <algorithm>
#include <string.h>
#include "objects.h"
#include <math.h>
#include <functional>
#include <numeric>
#include <iterator>

using namespace std;




bool EC::tipCellTest(void){
    
    //as the original tip cell test said there had to be 1.2x initmemagents no. which here is 740 which gives 148 left over to be involved in filopodia..
    //so as there would be slightly less spring agents than fil-memAgents in the previous model, as there are nodes in filopodia too we use 100
    //this remains the same even if cvell size is smaller with ECcross>1 as the no of filopdia agents should be the same..
    //if((VEGFR2tot>(float)VEGFR2norm/2.0f || VEGFR3tot>(float)VEGFR3norm/2.0f)&&(actinUsed>=257.0))return(true);
    if(((VEGFR2tot+VEGFR3tot)>((float)VEGFR2norm/2.0f+(float)VEGFR3norm/2.0f))&&(actinUsed>=257.0))return(true);
    //if((VEGFR2tot>(float)VEGFR2norm/2.0f)&&(actinUsed>=257.0))return(true);
    //LC - VEGFR2toVEGFR3// if((VEGFR3tot>(float)VEGFR3norm/2.0f)&&(VEGFR2tot>(float)VEGFR2norm/2.0f)&&(actinUsed>=257.0))return(true);
    else return(false);
    
}
//----------------------------------------------------------------------------------
/**
 * 
 * calculates total filopodia lengths on cell to see if overall maximum is reached
 * , if reached it will stop any further filopdia growth until some retract and more actin in available
 */
void EC::calcCurrentActinUsed(void){
   
    actinUsed=0;

    int i;
    for(i=0;i<nodeAgents.size();i++)
        if(nodeAgents[i]->FIL==TIP){
           
            actinUsed += nodeAgents[i]->FilLength(TIP);
        }
    
  
}
//------------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------------------------------------------------------------
int EC::calcVonNeighs(void){
    int i;
    int count=0;
    for(i=0;i<(int)nodeAgents.size();i++){
        nodeAgents[i]->vonNeu=false;
        if(nodeAgents[i]->checkNeighsVonForEnv()==true){
            count++;
            nodeAgents[i]->vonNeu=true;
            //nodeAgents[i]->labelled = true;
        }
    }
    
    for(i=0;i<(int)surfaceAgents.size();i++){
        surfaceAgents[i]->vonNeu=false;
        if(surfaceAgents[i]->checkNeighsVonForEnv()==true){
            count++;
            surfaceAgents[i]->vonNeu=true;
            //surfaceAgents[i]->labelled = true;
            
        }
    }
    VonNeighs = count+springAgents.size();
    return(count);
}

//----------------------------------------------------------------------------------------------------------------------------------------------
void EC::allocateProts(void){
    
    int i, j, div;
    int divJunction=0;
    int count=0;
    
    calcVonNeighs();
    div=VonNeighs;//nodeAgents.size()+springAgents.size()+surfaceAgents.size();

    

    //count how many agents are at the junction-----------------------------------------------
    for(i=0;i<(int)nodeAgents.size();i++){
        nodeAgents[i]->assessed=false;
        nodeAgents[i]->JunctionTest(true);
        if(nodeAgents[i]->junction==true) divJunction++;
        
    }
    for(i=0;i<(int)springAgents.size();i++){
        springAgents[i]->assessed=false;
        springAgents[i]->JunctionTest( true);
        if(springAgents[i]->junction==true) divJunction++;
    }
    for(i=0;i<(int)surfaceAgents.size();i++){
        surfaceAgents[i]->assessed=false;
        surfaceAgents[i]->JunctionTest( true);
        if(surfaceAgents[i]->junction==true) divJunction++;
        
    }
    
    
    /*int MfilTot = div-count;
     * int MnotFilTot = count;
     * float alpha = 1.3;
     * float left;
     * left = VEGFR2tot-(MfilTot*((float)(VEGFR2tot*alpha)/(float)div));*/
    //-----------------------------------------------------------------------------------------------------
    //set membrane prot levels for each cell----------------------------------------------------
    
    for(j=0;j<(int)nodeAgents.size();j++){
        
        
        
        
        if(nodeAgents[j]->vonNeu==true){
            nodeAgents[j]->VEGFR2=(float)VEGFR2tot/(float)div;
            nodeAgents[j]->VEGFR3=(float)VEGFR3tot/(float)div; //LC//
        }
        
        if(nodeAgents[j]->junction==true){
            
            nodeAgents[j]->Notch1=(float)NotchNorm/(float)divJunction;
            nodeAgents[j]->Dll4=(float)Dll4tot/(float)divJunction;
            
        }
        else{
            nodeAgents[j]->Notch1=0.0f;
            nodeAgents[j]->Dll4=0.0f;
        }
    }
    for(j=0;j<(int)springAgents.size();j++){
        
        
        
        
        springAgents[j]->VEGFR2=(float)VEGFR2tot/(float)div;
        springAgents[j]->VEGFR3=(float)VEGFR3tot/(float)div; //LC//
        
        //clustered VR-2 to filopodia
        //if(springAgents[j-nodeAgents.size()]->FIL!=NONE) springAgents[j-nodeAgents.size()]->VEGFR2=(float)(VEGFR2tot*alpha)/(float)div;
        //else springAgents[j-nodeAgents.size()]->VEGFR2=(float)left/(float)MnotFilTot;
        
        if(springAgents[j]->junction==true){
            springAgents[j]->Notch1=(float)NotchNorm/(float)divJunction;
            springAgents[j]->Dll4=(float)Dll4tot/(float)divJunction;
            
        }
        else{
            springAgents[j]->Notch1=0.0f;
            springAgents[j]->Dll4=0.0f;
        }
        
        
    }
    for(j=0;j<(int)surfaceAgents.size();j++){
        
        
        
        
        if(surfaceAgents[j]->vonNeu==true){
            surfaceAgents[j]->VEGFR2=(float)VEGFR2tot/(float)div;
            surfaceAgents[j]->VEGFR3=(float)VEGFR3tot/(float)div; //LC//
        }
        
        //clustered VR-2 to filopodia
        //if(springAgents[j-nodeAgents.size()]->FIL!=NONE) springAgents[j-nodeAgents.size()]->VEGFR2=(float)(VEGFR2tot*alpha)/(float)div;
        //else springAgents[j-nodeAgents.size()]->VEGFR2=(float)left/(float)MnotFilTot;
        
        if(surfaceAgents[j]->junction==true){
            surfaceAgents[j]->Notch1=(float)NotchNorm/(float)divJunction;
            surfaceAgents[j]->Dll4=(float)Dll4tot/(float)divJunction;
            
        }
        else{
            surfaceAgents[j]->Notch1=0.0f;
            surfaceAgents[j]->Dll4=0.0f;
        }
        
    }
    
    
    //--------------------------------------------------------------------------------------------------
}
//-------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------
void EC::NotchDelay(void){
    
    int i;
    vector<float>::iterator T=NotchDelayArray.begin();
    actNotCurrent=0.0f;
    
    //add new activeNotch value to delay array
    NotchDelayArray.push_back(activeNotchtot);
    
    //add the first element in the delay vector (has passed its delay time and can now have an effect)
    //to the end of the lasts vector - where it will have an effect until it leaves the lasts stack.
    NotchLastsArray.push_back(NotchDelayArray.front());
    
    //remove element from delay stack
    NotchDelayArray.erase(T);
    
    //remove frst element in lasts array, it nolonger has an effect on VEGFR2 levels.
    T=NotchLastsArray.begin();
    NotchLastsArray.erase(T);
    
    //update current amount of active notch that has an effect on VEGFR2 levels (the sum of the lasts array)
    
    for(i=0;i<actNot_VEGFR_lasts;i++){
        actNotCurrent=actNotCurrent+NotchLastsArray[i];
    }
    
    
}
//-------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------
void EC::VEGFRDelay(void){
    
    int i;
    //LC// actVEGFRcurrent=0.0f;
    actR2R2current=0.0f; //LC//
    actR2R3current=0.0f; //LC//
    actR3R3current=0.0f; //LC//
    //LC// vector<float>::iterator T=VEGFRDelayArray.begin();
    vector<float>::iterator T_R2R2=R2R2DelayArray.begin(); //LC//
    vector<float>::iterator T_R2R3=R2R3DelayArray.begin(); //LC//
    vector<float>::iterator T_R3R3=R3R3DelayArray.begin(); //LC//
    
    //LC// VEGFRDelayArray.push_back(activeVEGFRtot);
    R2R2DelayArray.push_back(activeR2R2tot); //LC//
    R2R3DelayArray.push_back(activeR2R3tot); //LC//
    R3R3DelayArray.push_back(activeR3R3tot); //LC//
    
    //add the first element in the delay vector (has passed its delay time and can now have an effect)
    //to the end of the lasts vector - where it will have an effect until it leaves the lasts stack.
    //LC// VEGFRlastsArray.push_back(VEGFRDelayArray.front());
    R2R2lastsArray.push_back(R2R2DelayArray.front()); //LC//
    R2R3lastsArray.push_back(R2R3DelayArray.front()); //LC//
    R3R3lastsArray.push_back(R3R3DelayArray.front()); //LC//
    
    //remove element from delay stack
    //LC// VEGFRDelayArray.erase(T);
    R2R2DelayArray.erase(T_R2R2); //LC//
    R2R3DelayArray.erase(T_R2R3); //LC//
    R3R3DelayArray.erase(T_R3R3); //LC//
    
    //remove frst element in lasts array, it nolonger has an effect on VEGFR2 levels.
    //LC// T=VEGFRlastsArray.begin();
    //LC// VEGFRlastsArray.erase(T);
    T_R2R2=R2R2lastsArray.begin(); //LC//
    R2R2lastsArray.erase(T_R2R2); //LC//
    T_R2R3=R2R3lastsArray.begin(); //LC//
    R2R3lastsArray.erase(T_R2R3); //LC//
    T_R3R3=R3R3lastsArray.begin(); //LC//
    R3R3lastsArray.erase(T_R3R3); //LC//
    
    //update current amount of active notch that has an effect on VEGFR2 levels (the sum of the lasts array)
    for(i=0;i<VEGFR_dll4_lasts;i++){
        //LC// actVEGFRcurrent=actVEGFRcurrent+VEGFRlastsArray[i];
        actR2R2current=actR2R2current+R2R2lastsArray[i]; //LC//
        actR2R3current=actR2R3current+R2R3lastsArray[i]; //LC//
        actR3R3current=actR3R3current+R3R3lastsArray[i]; //LC//
    }
}
//-------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------
void EC::GRN(void){
    
    //down-reg VEGFR2 via notch
    
    //LC - VEGFR2toVEGFR3// DLL4-VEGFRs affinity values
    affR2R2DLL4 = 1.0; //LC// 0.7-1.5
    affR2R3DLL4 = 0.67; //0.5; //LC// 0.1-0.8
    affR3R3DLL4 = 0.05; //LC// 0.0 ou 0.01
  
    VEGFR2tot=(VEGFR2norm)+NotchToVEGFR2*(actNotCurrent*sigma); //VEGFR2norm is now a EC specific param and scaled at config if mutant
    VEGFR3tot=(VEGFR3norm)+NotchToVEGFR3*(actNotCurrent*sigma); //VEGFR3norm is now a EC specific param and scaled at config if mutant

    if(VEGFR2tot<VEGFR2min) VEGFR2tot=VEGFR2min;
    if(VEGFR3tot<VEGFR3min) VEGFR3tot=VEGFR3min;

    //LC// Check how to combine the effects of the different dimers on DLL4
    actVEGFRcurrent = affR2R2DLL4 * actR2R2current + affR2R3DLL4 * actR2R3current + affR3R3DLL4 * actR3R3current; //LC//
    if(ANALYSIS_HYSTERESIS==true){
        if((this!=worldP->ECagents[0])&&(this!=worldP->ECagents[ECELLS-1]))
            Dll4tot+=(actVEGFRcurrent*delta);
        else{
        //for the fixed cells at the ends of the row, fix Dll4tot at current increment.
        Dll4tot = worldP->ECagents[1]->hyst.Current_Dll4_incremented_level;
        }
    }
    else Dll4tot+=(actVEGFRcurrent*delta);
    
    if(Dll4tot>MAX_dll4){
        Dll4tot=MAX_dll4;
    }
    
}
//---------------------------------------------------------------------------------

//---------------------------------------------------------------------------------
//get new protein totals from Magents
void EC::updateProteinTotals(void){
    
    int m;
    int uptoN=nodeAgents.size();
    int uptoS=springAgents.size();
    int uptoSu=surfaceAgents.size();
    
    //VEGFR2tot=0.0f;
    Dll4tot=0.0f;
    activeNotchtot=0.0f;
    //LC// activeVEGFRtot=0.0f;
    activeR2R2tot=0.0f;
    activeR2R3tot=0.0f;
    activeR3R3tot=0.0f;
    int junctionAgents=0;

    for(m=0;m<uptoN;m++){
        
        //LC-R2R2// activeVEGFRtot=activeVEGFRtot+nodeAgents[m]->VEGFRactive;
        activeR2R2tot=activeR2R2tot+nodeAgents[m]->R2R2active; //LC-R2R2//
        activeR2R3tot=activeR2R3tot+nodeAgents[m]->R2R3active; //LC-R2R3//
        activeR3R3tot=activeR3R3tot+nodeAgents[m]->R3R3active; //LC-R3R3//
        //if(nodeAgents[m]->junction==true){
        activeNotchtot=activeNotchtot+nodeAgents[m]->activeNotch;
        
        Dll4tot=Dll4tot+nodeAgents[m]->Dll4;
        //junctionAgents++;
        //}
    }
    for(m=0;m<uptoS;m++){
        
        //LC-R2R2// activeVEGFRtot=activeVEGFRtot+springAgents[m]->VEGFRactive;
        activeR2R2tot=activeR2R2tot+springAgents[m]->R2R2active; //LC-R2R2//
        activeR2R3tot=activeR2R3tot+springAgents[m]->R2R3active; //LC-R2R3//
        activeR3R3tot=activeR3R3tot+springAgents[m]->R3R3active; //LC-R3R3//
        //activeNotchtot=activeNotchtot+springAgents[m]->activeNotch;
        
        //Dll4tot=Dll4tot+springAgents[m]->Dll4;
    }
    for(m=0;m<uptoSu;m++){
        
        //LC-R2R2// activeVEGFRtot=activeVEGFRtot+surfaceAgents[m]->VEGFRactive;
        activeR2R2tot=activeR2R2tot+surfaceAgents[m]->R2R2active; //LC-R2R2//
        activeR2R3tot=activeR2R3tot+surfaceAgents[m]->R2R3active; //LC-R2R3//
        activeR3R3tot=activeR3R3tot+surfaceAgents[m]->R3R3active; //LC-R3R3//
        //if(surfaceAgents[m]->junction==true){
        activeNotchtot=activeNotchtot+surfaceAgents[m]->activeNotch;

        Dll4tot=Dll4tot+surfaceAgents[m]->Dll4;
        //junctionAgents++;
        //}
    }
    
    //worldP->dataFile<<Dll4tot<<" "<<activeNotchtot<<" "<<junctionAgents<<endl;
    NotchDelay();
    VEGFRDelay();
    //actNotCurrent = activeNotchtot;
    //actVEGFRcurrent = activeVEGFRtot;
}
//---------------------------------------------------------------------------------
/*void World::testNotchDll4(void){

int i,j,k,m;
int uptoN;
    int uptoS;
    int uptoSu;

    uptoN=nodeAgents.size();
    uptoS=springAgents.size();
    uptoSu=surfaceAgents.size();

    //VEGFR2tot=0.0f;
    Dll4tot=0.0f;
    activeNotchtot=0.0f;
    activeVEGFRtot=0.0f;
    int junctionAgents=0;

    int cells[4][4] = {0};

    for(i=0;i<ECagents.size();i++){

        for(m=0;m<uptoN;m++){
        if(nodeAgents[m]->junction==true){
            //need to check if it has aneigh of which cell
            if()
        }
        }
        for(m=0;m<uptoN;m++){
        if(surfaceAgents[m]->junction==true){

        }
        }
    }


        //LC-R2R2// activeVEGFRtot=activeVEGFRtot+nodeAgents[m]->VEGFRactive;
        activeVEGFRtot=activeVEGFRtot+nodeAgents[m]->R2R2active; //LC-R2R2//
        if(nodeAgents[m]->junction==true){
        activeNotchtot=activeNotchtot+nodeAgents[m]->activeNotch;

        Dll4tot=Dll4tot+nodeAgents[m]->Dll4;
        junctionAgents++;
        }
    }
    for(m=0;m<uptoS;m++){

        //LC-R2R2// activeVEGFRtot=activeVEGFRtot+springAgents[m]->VEGFRactive;
        activeVEGFRtot=activeVEGFRtot+springAgents[m]->R2R2active; //LC-R2R2//
        //activeNotchtot=activeNotchtot+springAgents[m]->activeNotch;

        //Dll4tot=Dll4tot+springAgents[m]->Dll4;
    }
    for(m=0;m<uptoSu;m++){

        //LC-R2R2// activeVEGFRtot=activeVEGFRtot+surfaceAgents[m]->VEGFRactive;
        activeVEGFRtot=activeVEGFRtot+surfaceAgents[m]->R2R2active; //LC-R2R2//
        if(surfaceAgents[m]->junction==true){
        activeNotchtot=activeNotchtot+surfaceAgents[m]->activeNotch;

        Dll4tot=Dll4tot+surfaceAgents[m]->Dll4;
        junctionAgents++;
        }
    }

    worldP->dataFile<<Dll4tot<<" "<<activeNotchtot<<" "<<junctionAgents<<endl;

    }*/


//------------------------------------------------------------------------------------------------
void EC::characterizeActNotchBoundaries(int which, int other){
    
    int i, j, y, n, k, m, p, s, x, t, B, flag;
   
    int sum=0;
    bool allow;
    int sumJunction=0;
    MemAgent* np;
    for(i=0;i<worldP->cellNeighbourhoods[which].size();i++){
        junctionSizes.push_back(0);
    }
    vector <EC*> neighbours;
    //hypothesis - why equal boundaries no offset fails more than equal boundaries offset - offset one is fighting on more fronts - or bigger fronts so has distributed its dll4 thinner all over as it were - so is essentially using less to inhibit the ones more likely to become tips - facing the vegf. We know that turing down dll4 helps the pattern to stabilise in high vegf...
    
    //group agents into which boundary they are in, with which cell,
    for(s=0;s<nodeAgents.size();s++){
        // - show how many neighbour cells it has  and how big each bounday is.. in percentage of agents..
        //and show total level of active notch on each boundary.
        if(nodeAgents[s]->junction==true){
             sumJunction++;
            np = nodeAgents[s];
            if((np->FIL==NONE)||(np->FIL==BASE)){
                neighbours.clear();
               
                i=np->Mx; j=np->My; k=np->Mz;
                
                //same layer
                for(x=0;x<26;x++){
                    
                    if(x==0){ m=i+1;n=j-1;p=k;}
                    else if(x==1){ m=i+1;n=j;p=k;  }
                    else if(x==2){ m=i+1;n=j+1;p=k;  }
                    else if(x==3){ m=i;n=j-1;p=k;  }
                    else if(x==4){ m=i;n=j+1;p=k;  }
                    else if(x==5){ m=i-1;n=j-1;p=k; }
                    else if(x==6){ m=i-1;n=j;p=k;  }
                    else if(x==7){ m=i-1;n=j+1;p=k;}
                    
                    //layer below
                    else if(x==8){m=i+1;n=j-1;p=k-1;}
                    else if(x==9) {m=i+1;n=j;p=k-1;}
                    else if(x==10){m=i+1;n=j+1;p=k-1;}
                    else if(x==11) {m=i;n=j-1;p=k-1;}
                    else if(x==12) {m=i;n=j+1;p=k-1;}
                    else if(x==13) {m=i-1;n=j-1;p=k-1;}
                    else if(x==14) {m=i-1;n=j;p=k-1;}
                    else if(x==15) {m=i-1;n=j+1;p=k-1;}
                    else if(x==16) {m=i;n=j;p=k-1;}
                    
                    //layer above
                    else if(x==17) {m=i+1;n=j-1;p=k+1;}
                    else if(x==18) {m=i+1;n=j;p=k+1;}
                    else if(x==19) {m=i+1;n=j+1;p=k+1;}
                    else if(x==20) {m=i;n=j-1;p=k+1;}
                    else if(x==21) {m=i;n=j+1;p=k+1;}
                    else if(x==22) {m=i-1;n=j-1;p=k+1;}
                    else if(x==23) {m=i-1;n=j;p=k+1;}
                    else if(x==24) {m=i-1;n=j+1;p=k+1;}
                    else {m=i;n=j;p=k+1;}
                    
                    
                    //-------------------------------
                    //toroidal only X
                    if(m>=xMAX) m=0;
                    if(m<0) m=xMAX-1;
                    //-------------------------------
                    
                    if(worldP->insideWorld(m, n, p)==true){
                        if(worldP->grid[m][n][p].type==M){
                            for(y=0;y<(int)worldP->grid[m][n][p].Mids.size();y++){
                                flag=0;
                                if(worldP->grid[m][n][p].Mids[y]->Cell!=this){
                                    //-----------------------------------------------------------------------
                                    if(other==-1){
                                    for(t=0;t<worldP->cellNeighbourhoods[which].size();t++){
                                        if(flag==0){
                                            if(worldP->grid[m][n][p].Mids[y]->Cell==worldP->cellNeighbourhoods[which][t]){
                                                flag=1;
                                                allow=true;
                                                for(B=0;B<neighbours.size();B++){
                                                    if(neighbours[B]==worldP->grid[m][n][p].Mids[y]->Cell) allow=false;
                                                }
                                                if(allow==true){
                                                    junctionSizes[t]++ ;
                                                    neighbours.push_back(worldP->grid[m][n][p].Mids[y]->Cell);
                                                }
                                                
                                            }
                                        }
                                    }
                                    if(flag==0) cout<<"new neighbour somehow?!"<<endl;
                                    }
                                    //----------------------------------------------------------------------
                                    //just find against a specific neighbour cell
                                    else{
                                        if(worldP->grid[m][n][p].Mids[y]->Cell==worldP->cellNeighbourhoods[which][other]){
                                                
                                                allow=true;
                                                for(B=0;B<neighbours.size();B++){
                                                    if(neighbours[B]==worldP->grid[m][n][p].Mids[y]->Cell) allow=false;
                                                }
                                                if(allow==true){
                                                    junctionSizes[t]++ ;
                                                    neighbours.push_back(worldP->grid[m][n][p].Mids[y]->Cell);
                                                }
                                                
                                            }
                                    }
                                    
                                }
                                
                            }
                        }
                        
                    }
                    
                }
            }
        }
    }
    
    
    for(i=0;i<worldP->cellNeighbourhoods[which].size();i++){
        cout<<junctionSizes[i]<<" ";
        sum+=junctionSizes[i];
    }
    cout<<endl;
    //cout<<"total: "<<sum<<"actual total(sum can be more tho as some agents are on two boundaries):"<<sumJunction<<endl;
    
}
//------------------------------------------------------------------------------------------------
