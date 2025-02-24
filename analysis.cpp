#include <math.h>
#include <iomanip>
#include "objects.h"
#include <fstream>
using namespace std;

//LC// What are these variables used for ??????
bool found_sp = false;
vector <int> current_order;
vector <int> shuffles;
vector <int> headShuffles_COM;
vector <int> orderFar;
vector <int> orderNear;
vector <int> stored_Far;
vector <int> stored_Near;
int oldHeadNode;
int overtakeTime = 1;
int frequencyHeads[ECELLS] = {0};
int frequencyHeadsDom[ECELLS] = {0};
int SnapShot_frequencyHeads[ECELLS] = {0};
int SnapShot_frequencyHeadsDom[ECELLS] = {0};
int battleTime = 0;
int domReighMutant=0;
int lastHead;
int totalTime = 0;
int snapShot_oldHead;
int dominant_reign = 0;
int countDomRei = 0;
int dom_head = -1;
int battling_reign = 0;
int battlers = 0;
int DOM_reigners[ECELLS] = {0};
int Battle_reigners[ECELLS] = {0};
int snapShot_count = 0;
bool snapShot_start = false;
int countActual = 0;
int TotalBattlingTime = 0;
int battles = 0;
int AILars = 0;
int countHET = 0;
int countHETact = 0;


//---------------------------------------------------------------------------------
//LC// Used in the .py files
//---------------------------------------------------------------------------------
std::vector< std::vector<float> > World::getGridSiteData()
{
    std::vector< std::vector<float> > retval;

//    for(int x = 0; x < xMAX; x++)
//    {
//        for(int y = 0; y < yMAX; y++)
//        {
//            for(int z = 0; z < zMAX; z++)
//            {
//                cout << "x: " << x << " y: " << y << " z: " << z << " type= " << grid[x][y][z].type << endl;
//            }
//        }
//    }

    for (int x = 0; x < xMAX; x++)
    {
        for (int y = 0; y < yMAX; y++)
        {
            for (int z = 0; z < zMAX; z++)
            {
                std::vector<float> gridSiteValues;
                gridSiteValues.push_back(x);
                gridSiteValues.push_back(y);
                gridSiteValues.push_back(z);

                if (grid[x][y][z].Eid != NULL)
                {
                    gridSiteValues.push_back(grid[x][y][z].Eid->VEGF);
                }
               else {
                    gridSiteValues.push_back(0);
                }

                float vegfTotal = 0;
                //LC// float vegfrTotal = 0; 
                float vegfr2Total = 0; //LC//
                float vegfr3Total = 0; //LC//
                //LC// float activeVegfrTotal = 0;
                float activeR2r2Total = 0;
                float activeR2r3Total = 0;
                float activeR3r3Total = 0;
                float dll4Total = 0;

                if (grid[x][y][z].type == E)
                {
                    gridSiteValues.push_back(0); //for environment type
                    gridSiteValues.push_back(grid[x][y][z].Fids.size());

                    for (int endothelialCellNumber = 0 ; endothelialCellNumber < ECELLS; endothelialCellNumber++)
                    {
                        vegfTotal = 0;
                        //LC// vegfrTotal = 0;
                        vegfr2Total = 0; //LC//
                        vegfr3Total = 0; //LC//
                        //LC// activeVegfrTotal = 0;
                        activeR2r2Total = 0; //LC//
                        activeR2r3Total = 0; //LC//
                        activeR3r3Total = 0; //LC//
                        dll4Total = 0;

                        if (grid[x][y][z].Fids.size() > 0) // Check if grid site contains filapodia
                        {
                            for (int filopodiaID = 0; filopodiaID < grid[x][y][z].Fids.size(); filopodiaID++)
                            {
                                // Check if filapodium belongs to this EC
                                if (grid[x][y][z].Fids[filopodiaID]->Cell == ECagents[endothelialCellNumber])
                                {
                                    vegfTotal += grid[x][y][z].Fids[filopodiaID]->SumVEGF;
                                    //LC// vegfrTotal += grid[x][y][z].Fids[filopodiaID]->VEGFR2;
                                    vegfr2Total += grid[x][y][z].Fids[filopodiaID]->VEGFR2; //LC//
                                    vegfr3Total += grid[x][y][z].Fids[filopodiaID]->VEGFR3; //LC//
                                    //activeVegfrTotal += grid[x][y][z].Fids[filopodiaID]->VEGFRactive;
                                    //activeVegfrTotal += grid[x][y][z].Fids[filopodiaID]->R2R2active; //LC-R2R2//
                                    activeR2r2Total += grid[x][y][z].Fids[filopodiaID]->R2R2active; //LC-R2R2//
                                    activeR2r3Total += grid[x][y][z].Fids[filopodiaID]->R2R3active; //LC-R2R3//
                                    activeR3r3Total += grid[x][y][z].Fids[filopodiaID]->R3R3active; //LC-R3R3//
                                    dll4Total += grid[x][y][z].Fids[filopodiaID]->Dll4;
                                }
                            }
                        }
                        gridSiteValues.push_back(vegfTotal);
                        //LC// gridSiteValues.push_back(vegfrTotal);
                        gridSiteValues.push_back(vegfr2Total); //LC//
                        gridSiteValues.push_back(vegfr3Total); //LC//
                        //LC// gridSiteValues.push_back(activeVegfrTotal);
                        gridSiteValues.push_back(activeR2r2Total); //LC//
                        gridSiteValues.push_back(activeR2r3Total); //LC//
                        gridSiteValues.push_back(activeR3r3Total); //LC//
                        gridSiteValues.push_back(dll4Total);
                    }
                }
                else if (grid[x][y][z].type == M)
                {
                    gridSiteValues.push_back(1); // for membrane type
                    gridSiteValues.push_back(grid[x][y][z].Fids.size() + grid[x][y][z].Mids.size());

                    for (int endothelialCellNumber = 0 ; endothelialCellNumber < ECELLS; endothelialCellNumber++)
                    {
                        vegfTotal = 0;
                        //LC// vegfrTotal = 0;
                        vegfr2Total = 0; //LC//
                        vegfr3Total = 0; //LC//
                        //LC// activeVegfrTotal = 0;
                        activeR2r2Total = 0; //LC//
                        activeR2r3Total = 0; //LC//
                        activeR3r3Total = 0; //LC//
                        dll4Total = 0;

                        if (grid[x][y][z].Fids.size() > 0)
                        {
                            for (int filopodiaID = 0; filopodiaID < grid[x][y][z].Fids.size(); filopodiaID++)
                            {
                                if (grid[x][y][z].Fids[filopodiaID]->Cell == ECagents[endothelialCellNumber])
                                {
                                    vegfTotal += grid[x][y][z].Fids[filopodiaID]->SumVEGF;
                                    //LC// vegfrTotal += grid[x][y][z].Fids[filopodiaID]->VEGFR2;
                                    vegfr2Total += grid[x][y][z].Fids[filopodiaID]->VEGFR2; //LC//
                                    vegfr3Total += grid[x][y][z].Fids[filopodiaID]->VEGFR3; //LC//
                                    //activeVegfrTotal += grid[x][y][z].Fids[filopodiaID]->VEGFRactive;
                                    //activeVegfrTotal += grid[x][y][z].Fids[filopodiaID]->R2R2active; //LC-R2R2//
                                    activeR2r2Total += grid[x][y][z].Fids[filopodiaID]->R2R2active; //LC-R2R2//
                                    activeR2r3Total += grid[x][y][z].Fids[filopodiaID]->R2R3active; //LC-R2R3//
                                    activeR3r3Total += grid[x][y][z].Fids[filopodiaID]->R3R3active; //LC-R3R3//
                                    dll4Total += grid[x][y][z].Fids[filopodiaID]->Dll4;
                                }
                            }
                        }
                        if (grid[x][y][z].Mids.size() > 0)
                        {
                            for (int memAgentID = 0; memAgentID < grid[x][y][z].Mids.size(); memAgentID++)
                            {
                                if (grid[x][y][z].Mids[memAgentID]->Cell == ECagents[endothelialCellNumber])
                                {
                                    vegfTotal += grid[x][y][z].Mids[memAgentID]->SumVEGF;
                                    //LC// vegfrTotal += grid[x][y][z].Fids[memAgentID]->VEGFR2;
                                    vegfr2Total += grid[x][y][z].Mids[memAgentID]->VEGFR2; //LC//
                                    vegfr3Total += grid[x][y][z].Mids[memAgentID]->VEGFR3; //LC//
                                    //activeVegfrTotal += grid[x][y][z].Mids[memAgentID]->VEGFRactive;
                                    //activeVegfrTotal += grid[x][y][z].Mids[memAgentID]->R2R2active; //LC-R2R2//
                                    activeR2r2Total += grid[x][y][z].Mids[memAgentID]->R2R2active; //LC-R2R2//
                                    activeR2r3Total += grid[x][y][z].Mids[memAgentID]->R2R3active; //LC-R2R3//
                                    activeR3r3Total += grid[x][y][z].Mids[memAgentID]->R3R3active; //LC-R3R3//
                                    dll4Total += grid[x][y][z].Mids[memAgentID]->Dll4;
                                }
                            }
                        }
                        gridSiteValues.push_back(vegfTotal);
                        //LC// gridSiteValues.push_back(vegfrTotal);
                        gridSiteValues.push_back(vegfr2Total); //LC//
                        gridSiteValues.push_back(vegfr3Total); //LC//
                        //LC// gridSiteValues.push_back(activeVegfrTotal);
                        gridSiteValues.push_back(activeR2r2Total); //LC//
                        gridSiteValues.push_back(activeR2r3Total); //LC//
                        gridSiteValues.push_back(activeR3r3Total); //LC//
                        gridSiteValues.push_back(dll4Total);
                    }
                }
                retval.push_back(gridSiteValues);
            }
        }
    }
    return retval;
}

//---------------------------------------------------------------------------------
//LC// Used in the .py files
//---------------------------------------------------------------------------------
std::vector< std::vector< std::vector < std::array<int,2> > > > World::getGridMapOfFilopodiaMovement()
 {
    vector < vector < vector < array<int, 2> > > > retval(xMAX, vector < vector < array<int, 2> > >(yMAX, vector < array<int, 2> >(zMAX, {0,0})));

    int totalExtenstions = 0;
    int totalRetractions = 0;

    for (EC* ec : ECagents)
    {
        for (std::array<int, 3> filExtension : ec->filopodiaExtensions)
        {
            totalExtenstions++;
            retval[filExtension[0]][filExtension[1]][filExtension[2]][0] += 1;
        }
        for (std::array<int, 3> filRetractions : ec->filopodiaRetractions)
        {
            totalRetractions++;
            retval[filRetractions[0]][filRetractions[1]][filRetractions[2]][1] += 1;
        }
    }

    return retval;
}

//---------------------------------------------------------------------------------
//LC// Used in the .py files
//---------------------------------------------------------------------------------
std::vector< std::vector< std::vector<float> > > World::getFilopodiaBaseLocationsAndTotalVegfr()
{
    std::vector< std::vector< std::vector<float> > > retval;

    std::vector< std::vector<float> > cell1Values;
    std::vector< std::vector<float> > cell2Values;

    if (this->timeStep == 800)
    {
        cout << "stop";
    }

    for (Filopodia* filopodia : this->filopodia)
    {
        if (!filopodia->retracted)
        {
            MemAgent* currentMemAgentInFilopodia = filopodia->base;
            if (currentMemAgentInFilopodia != nullptr)
            {
                float x = currentMemAgentInFilopodia->Mx;
                float y = currentMemAgentInFilopodia->My;
                float z = currentMemAgentInFilopodia->Mz;
                //float totalVegfr = currentMemAgentInFilopodia->VEGFRactive;
                float totalR2r2 = currentMemAgentInFilopodia->R2R2active; //LC-R2R2//
                float totalR2r3 = currentMemAgentInFilopodia->R2R3active; //LC-R2R3//
                float totalR3r3 = currentMemAgentInFilopodia->R3R3active; //LC-R3R3//

                //totalVegfr += currentMemAgentInFilopodia->VEGFRactive;
                totalR2r2 += currentMemAgentInFilopodia->R2R2active; //LC-R2R2//
                totalR2r3 += currentMemAgentInFilopodia->R2R3active; //LC-R2R3//
                totalR3r3 += currentMemAgentInFilopodia->R3R3active; //LC-R3R3//
                while (currentMemAgentInFilopodia->plusSite != nullptr)
                {
                    currentMemAgentInFilopodia = currentMemAgentInFilopodia->plusSite;
                    //totalVegfr += currentMemAgentInFilopodia->VEGFRactive;
                    totalR2r2 += currentMemAgentInFilopodia->R2R2active; //LC-R2R2//
                    totalR2r3 += currentMemAgentInFilopodia->R2R3active; //LC-R2R3//
                    totalR3r3 += currentMemAgentInFilopodia->R3R3active; //LC-R3R3//
                }

                std::vector<float> values;
                values.push_back(x); values.push_back(y); values.push_back(z); 
                values.push_back(totalR2r2); values.push_back(totalR2r3); values.push_back(totalR3r3);

                if (currentMemAgentInFilopodia->Cell == ECagents[0])
                {
                    cell1Values.push_back(values);
                }
                else if (currentMemAgentInFilopodia->Cell == ECagents[1])
                {
                    cell2Values.push_back(values);
                }
            }
            else {
                cout << "filopodia base is null what is going on ?????";
            }
        }
    }
    retval.push_back(cell1Values); retval.push_back(cell2Values);
    return retval;
}

void World::simulateTimestep(std::vector< std::vector<float> > cellIncrements)
{
    for (unsigned i = 0; i < cellIncrements.size(); i++)
    {
        EC* currentCell = ECagents[i];
        currentCell->VEGFR2tot += cellIncrements[i][0];
        currentCell->Dll4tot += cellIncrements[i][1];
        currentCell->Notchtot += cellIncrements[i][2];
        currentCell->VEGFR3tot += cellIncrements[i][3]; //LC// add a lign/column to cellIncrements to account for VEGFR3
    }
    simulateTimestep();
}


//---------------------------------------------------------------------------------
//LC// Used in the hysteresis simulation - should not concern us
//---------------------------------------------------------------------------------
bool Hysteresis::evaluate_hysteresis(ofstream& fileTo) {

    
    int i;
    
    //calc stability of current cell properties//
    
    //Dll4
    if((Cell->Dll4tot<stableDll4+Dll4_SigRange)&&(Cell->Dll4tot>stableDll4-Dll4_SigRange)){
        //stable, so increment timer
        //cout<<"stable.."<<stabilityTimer_latest<<endl;
        stabilityTimer_latest++;
        stabilityTimer_overall++;
    }
    else{
        //unstable so use current level as new start point to measure stability from and reset timer
        // cout<<"UNstable.."<<endl;
        stableDll4 = Cell->Dll4tot;
        stabilityTimer_latest=0;
    }
    //filopodia (use actin params)
        
    if((stabilityTimer_latest>=CELL_STABLE)||(stabilityTimer_overall>=bio_time_window)){
    //cell is deemed stable, so output results 
        
        cout<<"Stabilised! "<<Cell->Dll4tot<<"\t"<<stabilityTimer_overall<<endl;
        //cout<<"stable! "<<Current_Dll4_incremented_level<<endl;
        storeDll4.push_back(Cell->Dll4tot);
        storeTimes.push_back(stabilityTimer_overall);
        
        
        //increment fixed neigbour dll4 levels
        if(Current_Dll4_incremented_level==HYST_INCREMENT_MAX) direction=false;
        
        if(direction==false)//MAX_dll4)
                Current_Dll4_incremented_level-=HYST_INCREMENT;
        else 
            Current_Dll4_incremented_level+=HYST_INCREMENT;
        
        cout<<"hysteresis Dext = "<<Current_Dll4_incremented_level<<endl;
       
        
        //clear timers for next round.
        stabilityTimer_latest=0;
        stabilityTimer_overall=0;
        
    }
         if(Current_Dll4_incremented_level<0){ 
      
             //write data to file
             for(i=0;i<HYST_INCREMENT_MAX/HYST_INCREMENT;i++){
                 
             fileTo<<(i*HYST_INCREMENT)<<"\t"<<storeDll4[i]<<"\t"<<storeTimes[i]<<"\t"<<storeDll4[storeDll4.size()-1-i]<<"\t"<<storeTimes[storeDll4.size()-1-i]<<endl;
             }
             
             fileTo<<endl;
             fileTo<<storeDll4[HYST_INCREMENT_MAX/HYST_INCREMENT]<<"\t"<<storeTimes[HYST_INCREMENT_MAX/HYST_INCREMENT]<<endl;
             
             
             direction=true;
             
             storeDll4.clear();
             storeTimes.clear();
             
             return(false);
         }
         else return(true);

}


//---------------------------------------------------------------------------------
//LC// Print some values in a file - may be useful in the future but not necessary as it is
//---------------------------------------------------------------------------------
void World::printScores(ofstream& fileTo) {

    int i, j, k;
    int sumNewJunction = 0;
    int name, name2;
    int adjacentLoops = 0;
    int nonAdjacentLoops = 0;
    int c, s, uptoS;
    int uptoC = ECELLS;
    int count = 0;
    int flag;
    int neighCell = 0;
    EC* ecp;
    int Xsum = 0;
    int Patt = 1;
    int tipCount = 0;


    for (c = 0; c < uptoC; c++) {

        fileTo << ECagents[c]->VEGFR2tot << "\t";
    }

    fileTo << endl;
    /*	count=0;
            flag=0;
            ecp=ECagents[c];
            uptoS=ecp->stableArray.size();
		
		
            //------------------------------------------------------------------
            for(s=uptoS-1;s>=0;s--){
			
                    if((flag==0)&&(ecp->stableArray[s]==0)) flag=1;
			
                    if(flag==0) count++;
            }
		
            //cout<<"C= "<<count;
            //------------------------------------------------------------------
		
            //prints out how long it took to stabilize fully
            Xsum=Xsum+count;
		
		
            //print out how many tip cells and whether any are next to each other.
            //define tip cell as having over half full amount of VEGFR2 and more than initial amount of memAgents
		
            //fileTo<<ecp->VEGFR2tot<<"\t"<<ecp->actVEGFRArray[0];//<<ecp->Dll4tot<<"\t"<<ecp->actNotCurrent<<"\t";
		
            if(ecp->tipCellTest()==true){
                    //it is a tip cell
                    //fileTo<<"1\t";
                    //if its next to a tip cell score a bad point - not the pattern we want.
                    if(neighCell==1) Patt=0;
                    if(c==uptoC-1){
                            if(ECagents[0]->tipCellTest()==true) Patt=0;
    }
                    neighCell=1;
                    tipCount++;
            }
            else{
                    neighCell=0;
                    //fileTo<<"0\t";
            }
		
    }
    //fileTo<<endl;
    if(tipCount==0){Patt=0;}
	
    if(Patt==1){
            if((tipCount==4)||(tipCount==5)) Patt=1;
            else Patt=0;
    }
	
    fileTo<<100*(((float)Xsum/(float)uptoC)/(float)MAXtime)<<"\t"<<Patt<<"\t"<<tipCount<<endl;
     */



}


//---------------------------------------------------------------------------------
//LC// Evaluate the sand and pepper pattern to stop simulation
//---------------------------------------------------------------------------------

// Parameters init - SHOULD STAY OUT OF THE FUNCTION (init once and not each time locally in the function)
int patternHistory = 0;
bool patterned = false;

void World::evaluateSandP()
{
    // Init local variables
    int Patt = 0;
    int c, s, r;
    int count;
    int flag = 0;
    int neighCell = 0;

    int uptoS = ECELLS;
    int tipCount = 0;
    int cont = 0;
    EC* ecp;

    // Loop over all EC agents
    for (c = 0; c < ECELLS; c++)
    {
        // Collect current EC agent
        ecp = ECagents[c];

        if (ecp->tipCellTest())
        {
            // Increase tip counter if cell is tip
            tipCount++;
            // Break if a tip cell has a neighbour also tip - check tipCount if not
            if (neighCell == 1)
            {
                Patt = 0;
                //break;
            }
            // Condition for pattern: 4 or 5 tip cells
            else if (tipCount == 4 || tipCount == 5)
            {
                Patt = 1;
            }
            neighCell = 1;
        }
        else 
        {
            neighCell = 0;
        }
    }

    tipCountSoso = tipCount;

    // WRITE TIP COUNT IN FILE
    // TIP_evolution_file << tipCount << endl;
    if (timeStep != MAXtime+1)
    {
        RUNSfile << tipCount << endl;
    }

    if (Patt == 1)
        patternHistory++;
    else
        patternHistory = 0;

    if (patternHistory == 150) {
        cout << "patterned!" << endl;
        patterned = true;
        RUNSfile << timeStep << endl;
        TIMETOPATTERNfile << timeStep << endl;
        TIPAMOUNTfile << tipCount << endl;
        timeStepPattern = timeStep;
        timeStep = MAXtime;
    }

    if (timeStep == MAXtime - 1 && !patterned)
    {
        RUNSfile << "-1" << endl;
        TIMETOPATTERNfile << "-1" << endl;
        TIPAMOUNTfile << tipCount << endl;
    }
}


//---------------------------------------------------------------------------------
/**
 * Output the current total membrane grid sites and total VEGF level in world
 **/
//LC// Used in timestep creation only for now - may be useful if other VEGF are introduced
//---------------------------------------------------------------------------------
void World::calcEnvVEGFlevel(void)
{

    int i, j, k;
    float sum = 0.0f;

    int count = 0; // env agents (??)
    int countM = 0; // "gridded" mem agents (??)

    for (i = 0; i < xMAX; i++)
    {
        for (j = 0; j < yMAX; j++)
        {
            for (k = 0; k < zMAX; k++)
            {
                if ((grid[i][j][k].type == E) && (grid[i][j][k].Eid->VEGF > 0.0f))
                {
                    sum += grid[i][j][k].Eid->VEGF;
                    count++;
                }
                if ((grid[i][j][k].type == M) && (grid[i][j][k].Mids.size() >= 1))
                {
                    countM++;
                }
            }
        }
    }

    cout << "VEGF total " << sum << "\t" <<count<<"gridded membrane" << countM << endl;

}


//---------------------------------------------------------------------------------
//LC// Output protein levels -> may be adapted and used in our analyses
//---------------------------------------------------------------------------------
void World::output_cell_protlevels(ofstream& fileTo)
{
    // Init local variables
    int i,j;
    int filTokSum = 0;
    EC* ecp;
    
    // Loop over all EC agents
    for(i=0;i<ECagents.size();i++)
    {
        // Reset filTokSum and collect current EC agent
        filTokSum = 0;
        ecp = ECagents[i];

        // Tally up filTokens
        for(j=0;j<ecp->nodeAgents.size();j++)
        {
            filTokSum += ecp->nodeAgents[j]->filTokens;
        }

        // Output initial protein values at first timestep - output current proteins levels and activation otherwise
        if(timeStep==0)
        {
            cout<<"CEll: "<<i<<" VEGFR2norm: "<<ecp->VEGFR2norm<<" VEGFR3norm: "<<ecp->VEGFR3norm<<" Vsink: "<<ecp->Vsink<<endl;
        }
        else
        {
            cout<<timeStep<<" CEll: "<<i<<" actinUsed "<<ecp->actinUsed<<" fiToks "<<filTokSum<<" Dll4tot:"<<ecp->Dll4tot<<" activeR2R2 "<<ecp->activeR2R2tot<<" activeR2R3 "<<ecp->activeR2R3tot<<" activeR3R23 "<<ecp->activeR3R3tot<<" VEGFR2norm: "<<ecp->VEGFR2norm<<" VEGFR3norm: "<<ecp->VEGFR3norm<<" Vsink: "<<ecp->Vsink<<endl; 
        }

        // Collect protein levels in file if necessary (not used for now)
        //fileTo<<ecp->Dll4tot<<"\t"<<ecp->activeVEGFRtot<<"\t";
    }
    
    // End line in console or file
    cout<<endl; //fileTo<<endl;
}


//---------------------------------------------------------------------------------
//LC// Used to track cell motion - should not be used in our case
//---------------------------------------------------------------------------------
void EC::store_cell_COM(void)
{
//centre of mass of cell position - though not toroidal
    float aveX=0;
    float aveY=0;
    float aveZ=0;
    int i;
    int count=0;
    Coordinates* COM = new Coordinates();

    int uptoN=nodeAgents.size();

    for(i=0;i<uptoN;i++)
    {
        if(nodeAgents[i]->minusSite==NULL)
        {
        aveX+=nodeAgents[i]->Mx;
        aveY+=nodeAgents[i]->My;
        aveZ+=nodeAgents[i]->Mz;
        count++;
        }
    }

    COM->x = aveX/(float)count;
    COM->y = aveY/(float)count;
    COM->z = aveZ/(float)count;

    COMstore.push_back(COM);
}


//---------------------------------------------------------------------------------
//LC// Should not be used in our case, ANALYSIS_JTB_SP_PATTERN==false in objects.h
//---------------------------------------------------------------------------------
void EC::calcStability(void)
{
    // if its significantly different from current stable amount, score 0 and change stableVEGFR2 to current VEGFR2
    if( (VEGFR2tot>stableVEGFR2+SigRange) || (VEGFR2tot<stableVEGFR2-SigRange) )
    {
        stableArray.push_back(0);
        stableVEGFR2=VEGFR2tot;
    }
    else 
    {
        stableArray.push_back(1);
    } 
}