#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <iterator>
#include "objects.h"
#include <string.h>
#include <time.h>
#include <algorithm>
#include <math.h>
#include <chrono>


#if GRAPHICS
#include "display.h"

#ifdef __linux__

#include <GL/glut.h>

#include <GL/glu.h>

#include <GL/gl.h>

#elif __APPLE__

#include <GLUT/glut.h>

#include <OpenGL/glu.h>

#include <OpenGL/gl.h>

#endif
#endif

using namespace std;
//using std::random_shuffle;

// General variables
World* WORLDpointer;
ofstream RUNSfile;
ofstream TIMETOPATTERNfile;
ofstream TIPAMOUNTfile;
// ofstream TIP_evolution_file;
// ofstream DLL4_evolution_file;
// ofstream VEGFR2_evolution_file;
// ofstream VEGFR3_evolution_file;
// ofstream activeNotch_evolution_file;
// ofstream activeR2R2_evolution_file;
// ofstream activeR2R3_evolution_file;
// ofstream activeR3R3_evolution_file;
int memINIT;
char fname[200];
float actinMax = 512;

// Model alterations
float intersoso;

// GRN Signalling pathways variables
float delta = 2.0; //2.0f normal
float sigma = 15;// 10.35f; //15 normal JTB setup
float NotchNorm;
float MAX_dll4;
float VEGFR2NORM; //total of receptors it will maintain if all else is equal - divides out to M agents
float VEGFR3NORM; //LC//
float VEGFR3scale;
float VEGFR2min;
float VEGFR3min; //LC//
float DLL4initscale;

float NotchToVEGFR2;
float NotchToVEGFR3;

// Filopodia elongation - Distributions parameters
float beta_R2R3;
float beta_R3R3;

// HUVEC or LEC
int lymphatic;

// ENVIRONMENT SETUP
float VEGFconc = 0.8f; //for uniform VEGF above a vessel JTB 2008
float horV = 0.04; //for NCB 2010 blind ended sprout radiating gradient
float perpV = 0.04; //for NCB 2010 blind ended sprout radiating gradient
float HorCutOff = ECwidth*ECpack - (ECwidth / 2.0f); //for NCB 2010 blind ended sprout radiating

float VconcST = 0.04;
float VconcSTMACRO = 0.15f; // Macrophage point source for PLoS CB 2009

// Rearrangement CPM module
CPM_module* diffAd;
int MCS = 8000;
float M1_neta = 200.0f; //M1 differential adhesion neta parameter value (used to be called NUMERATOR and set to double this amount in previous code) to determine diffferential adhesion , or =0 for  all weakly adhesive, =5000 for all strongly adhesive 
float M2_lambda = 200.0f;

//remove! //LC// wtf?
bool MEM_LEAK_OCCURRING = false; //core removal

// Hysteresis related variables
bool continue_hysteresis;
float dll4_SIG = 7.0f;
float FIL_VARY = 2;
float EPSILON = 0.9;
float tokenStrength = 1;
float FILTIPMAX = 15;
int FIL_SPACING = 2;
int run_number = 1;
int GRADIENT = STEADY;
float randFilExtend = -1;
float RAND_FILRETRACT_CHANCE = -1;
long long seed = -1;
mt19937 g;
uniform_real_distribution<double> dist = uniform_real_distribution<double>(0, NEW_RAND_MAX);

// Paper specific parameters - Junctional offset simulations as in PLoS CB 2009
int Junct_arrange = UNEQUAL_NEIGHS;
float CellPosOffset;

//------------------------------------------------------------------------------
//LC// maybe should stay
//------------------------------------------------------------------------------
//#define BAHTI_ANALYSIS true
#if BAHTI_ANALYSIS
template <typename... Args>
using overload_cast_ = pybind11::detail::overload_cast_impl<Args...>;
PYBIND11_MODULE(springAgent, m) {
    py::class_<World>(m, "World")
            .def(py::init<float, float, int, float, float, float, int, float, float, long long>(), "World constructor: World(float epsilon, float vconcst, int gradientType, float filConstNorm, float filTipMax, float tokenstrength, int filSpacing, float randFilExtension, float randFilRetract, long long seed)")
            .def("simulateTimestep", overload_cast_<>()(&World::simulateTimestep), "Simulate one timestep in MemAgent-Spring Model")
            .def("simulateTimestep", overload_cast_<std::vector<std::vector<float>>>()(&World::simulateTimestep), "Simulate one timestep in MemAgent-Spring Model with 2d array of protein value changes for each cell.")
            .def_readwrite("timeStep", &World::timeStep)
            .def_readonly("gridXDimensions", &World::gridXDimensions)
            .def_readonly("gridYDimensions", &World::gridYDimensions)
            .def_readonly("gridZDimensions", &World::gridZDimensions)
            .def("getGridSiteData", &World::getGridSiteData, "Get grid site data in 2D array. Each index in first array represents a grid site then the second array contains X, Y, Z, grid site vegf, Grid Type (0 = environment, 1 = membrane), Number of sub agents (filopodia/memagents), Cell1 vegf (the amount of vegf the cell can see in neighbour grid sites from that grid site), Cell1 vegfr, Cell1 active vegfr, cell2 vegf, cell2 vegfr, cell2 active vegfr, cell3... etc etc")
            .def("getFilopodiaBaseLocationsAndTotalVegfr", &World::getFilopodiaBaseLocationsAndTotalVegfr, "returns a 3 dimensional vector. first dimension is cell 1 then cell 2. second dimension contains vectors for all the filopodia in the cell. 3rd dimension is the filopodia values, these are: 0. base x coordinate, 1. base y coordinate, 2. base z coordinate, 3. total active vegfr in the filopodia")
            .def("getGridMapOfFilopodiaMovement", &World::getGridMapOfFilopodiaMovement, "returns 3d vector grid reflecting grid in simulation. The inner most vector contains an array of size 2. The first value in the array represents how many filopodia extensions into that grid site have occured on the current timestep. Second value represents the amount of retractions.");
}
#endif
//---------------------------------------------------------------------------------
//LC// Read arguments
//---------------------------------------------------------------------------------
void readArgs(int argc, char * argv[])
{
    if (argc > 1)
    {
        run_number = atoi(argv[1]);
    }
    if (argc > 2)
    {
        EPSILON = atof(argv[2]);
        VconcST = atof(argv[3]);
        GRADIENT = atoi(argv[4]);
        FIL_VARY = atof(argv[5]);
        FILTIPMAX = atof(argv[6]);
        tokenStrength = atof(argv[7]);
		FIL_SPACING = atof(argv[8]);
        VType = atoi(argv[9]);
        intersoso = atof(argv[10]);
        beta_R2R3 = atof(argv[11]);
        beta_R3R3 = atof(argv[12]);
        VEGFR3scale = atof(argv[13]);
        DLL4initscale = atof(argv[14]);
        NotchToVEGFR2 = atof(argv[15]);
        NotchToVEGFR3 = atof(argv[16]);
        lymphatic = atoi(argv[17]);
        if (argc > 17)
        {
            randFilExtend = atof(argv[18]);
            if (randFilExtend >= 0 && randFilExtend <= 1)
            {
                EPSILON = 0;
            }
            RAND_FILRETRACT_CHANCE = atof(argv[19]);
            if (argc > 20)
            {
                seed = stoll(argv[20]);
            }
        }
        VEGFconc = VconcST;
    }
}


//---------------------------------------------------------------------------------
//LC// MAIN LOOP YOUHOU //LC//
//---------------------------------------------------------------------------------
int main(int argc, char * argv[])
{

    // Read arguments and create statistics file
	char statistics_file_buffer[500];
    readArgs(argc, argv);

    if (ANALYSIS_HYSTERESIS)
    {
		sprintf(statistics_file_buffer,
				"statistics_hysteresis_filvary_%g_epsilon_%g_VconcST%g_GRADIENT%i_FILTIPMAX%g_tokenStrength%g_FILSPACING%i_randFilExtend%g_randFilRetract%g_seed%lld_run%i_.csv",
				double(FIL_VARY), double(EPSILON), VconcST, GRADIENT, FILTIPMAX, tokenStrength, FIL_SPACING,
				randFilExtend, RAND_FILRETRACT_CHANCE, seed, run_number);
	} 
    else if (ANALYSIS_TIME_TO_PATTERN)
    {
		sprintf(statistics_file_buffer,
				"statistics_time_to_pattern_filvary_%g_epsilon_%g_VconcST%g_GRADIENT%i_FILTIPMAX%g_tokenStrength%g_FILSPACING%i_randFilExtend%g_randFilRetract%g_seed%lld_run%i_.csv",
				double(FIL_VARY), double(EPSILON), VconcST, GRADIENT, FILTIPMAX, tokenStrength, FIL_SPACING,
				randFilExtend, RAND_FILRETRACT_CHANCE, seed, run_number);
    }
    else
    {
		sprintf(statistics_file_buffer,
				"statistics_msm_filvary_%g_epsilon_%g_VconcST%g_GRADIENT%i_FILTIPMAX%g_tokenStrength%g_FILSPACING%i_randFilExtend%g_randFilRetract%g_seed%lld_run%i_.csv",
				double(FIL_VARY), double(EPSILON), VconcST, GRADIENT, FILTIPMAX, tokenStrength, FIL_SPACING,
				randFilExtend, RAND_FILRETRACT_CHANCE, seed, run_number);
    }

	create_statistics_file(statistics_file_buffer);

    // Initialise timestep, print in file and output simulation parameters in console
	std::time_t start_time = get_current_time();
	string start_time_string = format_time_string(start_time, true);
	write_to_statistics_file(statistics_file_buffer, start_time_string);

	cout << "Current time: " << std::ctime(&start_time);
    cout << "ECPACK: " << ECpack << endl;
    cout << "GRAPHICS: " << GRAPHICS << endl;
    cout << "bahti analysis: " << BAHTI_ANALYSIS << " @@ time to pattern analysis: " << ANALYSIS_TIME_TO_PATTERN << " @@  hysteresis analysis: " << ANALYSIS_HYSTERESIS << endl;
    cout << "ECELLS: " << ECELLS << endl;
    cout << "Epsilon: " << EPSILON << endl;
    cout << "VconcST: " << VconcST << endl;
    cout << "gradient: " << GRADIENT << endl;
    cout << "FIL_VARY (filconstnorm): " << FIL_VARY << endl;
    cout << "FILTIPMAX: " << FILTIPMAX << endl;
    cout << "tokenStrength: " << tokenStrength << endl;
	cout << "FIL_SPACING: " << FIL_SPACING << endl;
	cout << "randFilExtension: " << randFilExtend << endl;
	cout << "RAND_FILRETRACT_CHANCE: " << RAND_FILRETRACT_CHANCE<< endl;

    // Initialise output files names ans open files
    char outfilename[500];
    char tipamountfilename[500];
    char timetopatternfilename[500];
    // char tipamountevolution[500];

    if (ANALYSIS_HYSTERESIS) {
        cout << "running bistability analysis" << endl;
        sprintf(outfilename, "analysis_hysteresis_filvary_%g_epsilon_%g_VconcST%g_GRADIENT%i_FILTIPMAX%g_tokenStrength%g_FILSPACING%i_randFilExtend%g_randFilRetract%g_seed%lld_VType%i_run%i.txt", double(FIL_VARY), double(EPSILON), VconcST, GRADIENT, FILTIPMAX, tokenStrength, FIL_SPACING, randFilExtend, RAND_FILRETRACT_CHANCE, seed, VType, run_number);
    }
    else if (ANALYSIS_TIME_TO_PATTERN) {
        cout << "running time to pattern analysis" << endl;
        sprintf(outfilename, "time_to_pattern_filvary_%g_eps_%g_VconcST%g_GRAD%i_FILTIPMAX%g_tokenStrength%g_FILSPACING%i_randFilExt%g_randFilRetract%g_seed%lld_VType%i_intersoso%g_betaR2R3%g_betaR3R3%g_VEGFR3scale%g_DLL4initscale%g_NotchToR2%g_NotchToR3%g_lymphatic%i_run%i.txt", double(FIL_VARY), double(EPSILON), VconcST, GRADIENT, FILTIPMAX, tokenStrength, FIL_SPACING, randFilExtend, RAND_FILRETRACT_CHANCE, seed, VType, intersoso, beta_R2R3, beta_R3R3, VEGFR3scale, DLL4initscale, NotchToVEGFR2, NotchToVEGFR3, lymphatic, run_number);
    }
    else {
        cout << "analysis must either be ANALYSIS_HYSTERESIS or ANALYSIS_TIME_TO_PATTERN.. aborting run";
        sprintf(outfilename, "testforpybind");
        //exit(1);
    }

    sprintf(timetopatternfilename, "output/time_to_pattern_concentration%g_parameter%g_run%i.txt", VconcST, intersoso, run_number);
    sprintf(tipamountfilename, "output/number_of_tip_concentration%g_parameter%g_run%i.txt", VconcST, intersoso, run_number);
    // sprintf(tipamountfilename, "output/tip_evolution_concentration%g_parameter%g_run%i.txt", VconcST, intersoso, run_number);

    cout << "output file name: " << outfilename << endl;

    RUNSfile.open(outfilename);
    TIMETOPATTERNfile.open(timetopatternfilename, ios::app);
    TIPAMOUNTfile.open(tipamountfilename, ios::app);
    // TIP_evolution_file.open("output/tip-evolution.txt");
    // DLL4_evolution_file.open("output/dll-evolution.txt");
    // VEGFR2_evolution_file.open("output/vegfr2-evolution.txt");
    // VEGFR3_evolution_file.open("output/vegfr3-evolution.txt");
    // activeNotch_evolution_file.open("output/active-notch-evolution.txt");
    // activeR2R2_evolution_file.open("output/active-r2r2-evolution.txt");
    // activeR2R3_evolution_file.open("output/active-r2r3-evolution.txt");
    // activeR3R3_evolution_file.open("output/active-r3r3-evolution.txt");

    // Create main world variable
    World* world = new World();
    WORLDpointer = world;
	
    //LC// what is it used for?? //LC//
    for (int i=0; i<16; ++i)
            new_rand();
    std::vector<int> testNumber = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    new_random_shuffle(testNumber.begin(), testNumber.end());
    //LC// what is it used for?? //LC//

    // Create main CPM variable (may not be used in our simulations I think)
    diffAd = new CPM_module(world);

    // Start graphics GUI if asked
#if GRAPHICS
    //main display function - simulating the model is called within here
    displayGlui(&argc, argv);
    glutMainLoop();
#else

    // LAUNCH MAIN SIMULATION
    world->runSimulation();

    // Write in output files
    /// Get end time, and calculate elapsed time -> add these to results file.
	std::time_t end_time = get_current_time();
	std::cout << "End time: " << std::ctime(&end_time) << std::endl;

	string end_time_string = format_time_string(end_time, false);
	write_to_statistics_file(statistics_file_buffer, end_time_string);

	long elapsed_time = end_time - start_time;
	string elapsed_time_string = "Elapsed time, " + to_string(elapsed_time);
	write_to_statistics_file(statistics_file_buffer, elapsed_time_string);

#endif

    // Close output files
    RUNSfile.close();
    TIMETOPATTERNfile.close();
    TIPAMOUNTfile.close();
}


//------------------------------------------------------------------------------------------
//LC// Main function - launching the simulation
//------------------------------------------------------------------------------------------
void World::runSimulation()
{
    // Simulate timesteps until MAXtime
    while (timeStep <= MAXtime)
    {
        // Print timestep in console every 50 timesteps
        if (timeStep % 50 == 0)
        {
            cout << "timestep " << timeStep << ".. " << MAXtime - timeStep << " left" << endl;
        }

        // Simulate current timestep
        simulateTimestep();

        // Run simulation analysis
        if (ANALYSIS_HYSTERESIS)
        {
            hysteresisAnalysis();
        }
        else if (ANALYSIS_TIME_TO_PATTERN)
        {
            evaluateSandP();
        }

        // Memory issues
        if(MEM_LEAK_OCCURRING)
        {
            timeStep = MAXtime;
            RUNSfile<<"MEMORY LEAKED!!!...quit run"<< endl;
            cout << "MEMORY LEAKED!!!...quit run" << endl;
            MEM_LEAK_OCCURRING = false;
        }

        // Go twice to next line in output file if we reached MAXtime
        if (timeStep == MAXtime)
        {
            RUNSfile << endl << endl;
        }

        /* //LC// MAY BE USED IN THE FUTURE //LC//
        // Get gridsite data at specific timeStep
        if (timeStep ==3)
        {
            getGridSiteData();
        }
        printScores(RUNSfile, RUNSfile2, RUNSfile3);
        */ //LC// MAY BE USED IN THE FUTURE //LC//
    }
    cout << "end of run simulation" << endl;
}

//------------------------------------------------------------------------------------------
//LC// Step simulation function (called once)
// --> creationTimeStep()
// --> updateMemAgents()
// --> updateECagents()
//------------------------------------------------------------------------------------------
void World::simulateTimestep()
{
    int movie = 0; //LC// Seems useless (will always be 0...?)
    timeStep++;

    // Initialise everything at first timestep - simulate timestep otherwise
    //LC// I agree with the following comments but let's keep it like that for now
    //TODO: maybe move this out of simulate timestep? bit misleading that its in here
    //could just call creation timestep func from here.. and have timesteps start from zero instead of -1
    if (timeStep == 0)
    {
        cout << "Creation timestep... initialising everything" << endl;
        creationTimestep(movie);
    }
    else
    {
        // Clear extensions and retractions variables (??figure out why??)
        for (EC* ec : ECagents)
        {
            ec->filopodiaExtensions.clear();
            ec->filopodiaRetractions.clear();
        }

        // Update each mem agent
        updateMemAgents();

        // Run CPM if conditions met - //LC// should not happen in our simulations
        if ( (timeStep > TIME_DIFFAD_STARTS) && REARRANGEMENT)
        {
            diffAd->run_CPM();
        }

        // Update each ec agent
        updateECagents();

        //LC// MOVIE MAKING HERE ??
        //movieMaking(movie);
    }
}

//------------------------------------------------------------------------------------------
//LC// FIRST step simulation function (called once, at first timestep)
//------------------------------------------------------------------------------------------
void World::creationTimestep(int movie)
{
    // Create macrophages if flag true - //LC// should not happen in our case
    if (MACROS > 0)
    {
        createMacrophages();
    }

    // Create EC agents depending on setup - //LC// we are using CELL_SETUP==1 --> vessel
    /** create EC cells spring mesh and memAgents within continuous space **/
    if ((CELL_SETUP == 2) || (CELL_SETUP == 1))
    {
        //blind ended sprout (NCB and rearrangement)
        //vessel setup (JTB and PLoS)
        createECagents(Junct_arrange);
        connectMesh();
    }
    else if (CELL_SETUP == 3)
    {
        createMonolayer();
    }
    else if (CELL_SETUP == 4)
    {
        create_3D_round_cell();
    }

    // Place agents onto gridded lattice
    for (int j = 0; j < (int) ECagents.size(); j++)
    {
        ECagents[j]->gridAgents();
    }

    // Collect the initial number of mem agents and output the number in the console
    /** set the memInit value if needed for watching cell growth and tip cell quantification **/
    memINIT = ECagents[0]->nodeAgents.size()+ ECagents[0]->surfaceAgents.size();
    cout << "memInit" << memINIT << endl;

    // Create environment
    createEnvironment();
    cout <<"created environment"<<endl;

    //LC// Seems unused - confirm when testing CELL_SETUP==2 ?
    ///TODO: ask kate if this still needs to be in here?
    if (CELL_SETUP == 2) {
        /*if (BLINDENDED_SPROUT == true) {
            //sew up leading front if blindended sprout setup
            for (i = 0; i < ECagents[ECELLS - 1]->nodeAgents.size(); i++) {
                if (ECagents[ECELLS - 1]->nodeAgents[i]->labelledBlindended == true) {

                    ECagents[ECELLS - 1]->nodeAgents[i]->moveAgent(ECagents[ECELLS - 1]->nodeAgents[i]->Mx, ECagents[ECELLS - 1]->nodeAgents[i]->My, (float) (vesselCentreZ - (vesselRadius - 1)), true);

                }
            }
            for (i = 0; i < ECagents[ECELLS - 2]->nodeAgents.size(); i++) {
                if (ECagents[ECELLS - 2]->nodeAgents[i]->labelledBlindended == true) {

                    ECagents[ECELLS - 2]->nodeAgents[i]->moveAgent(ECagents[ECELLS - 2]->nodeAgents[i]->Mx, ECagents[ECELLS - 2]->nodeAgents[i]->My, (float) (vesselCentreZ - (vesselRadius - 1)), true);


                }
            }
        }*/
    }

    if (ASTRO == RETINA)
    {
        create_astro_retina_section();
    }
    else if (ASTRO != NONE)
    {
        createAstrocytes();
    }

    // Create blood inside the vessels (== label inside and set VEGF=0)
    createBlood(); //labels the interior of vessels but not otherwise involved

    // Set VEGF levels according to input parameters
    /*give Env objects correct VEGF level depending on chosen gradients*/
    setInitialVEGF();

    // Split EC protein levels into its mem agents
    /*divide out cells overall levels of proteins to their memAgents once created*/
    for (int j = 0; j < (int) ECagents.size(); j++)
    {
        ECagents[j]->allocateProts();
    }

    // Define mem agents which to receptor activations
    /*define exposed membrane agents as those with a face adjacent to env not vertex (von neu neighbours)
    only these are flagged to do receptor interactions (required to give matching behaviour when scaling grid)*/
    label_env_exposed_von_neu_agents();

    // Output total VEGF, and the number of env agents and gridded (?) mem agents (Note: flag always true)
    if (ANALYSIS_TOTALVEGF_TOTAL_MEMBRANE)
    {
        calcEnvVEGFlevel();
    }

    //LC// CPM module set up - should not be used in our case 
    ///on first timestep this sets up the CPM module
    if (REARRANGEMENT)
    {
        diffAd->run_CPM();
    }

    // Output EC initial protein levels in console - always true
    if (ANALYSIS_PROTLEVELS)
    {
        output_cell_protlevels(dataFile);
    }

    //LC// Still do not get this movie thin (:
    /**
     * TODO: sort this movie thing out to hopw it was before.. only used in graphics on, so maybe hav
     * separate setup func to this
    **/
    if (movie == 1) system(" mkdir /tmp/movie2; rm -vf /tmp/movie2/*");

    // Clean and (re-)create movie directory
    system("mkdir movie; rm -vf movie/*");
    cout << "Creation complete" << endl;
}

//------------------------------------------------------------------------------------------
//LC// Hysteresis analysis (not our focus for now)
//------------------------------------------------------------------------------------------
void World::hysteresisAnalysis()
{
    if (timeStep == 1)
    {
        ECagents[1]->hyst.stableDll4 = ECagents[1]->Dll4tot;
        ECagents[1]->hyst.stableActin = ECagents[1]->actinUsed;
    }
    continue_hysteresis = ECagents[1]->hyst.evaluate_hysteresis(RUNSfile);
    if (!continue_hysteresis)
    {
        timeStep = MAXtime+1;
        RUNSfile<<endl<<endl;
     }
}

//------------------------------------------------------------------------------------------
//LC// mem agents update (called at each timestep)
//------------------------------------------------------------------------------------------
/**
 * Asynchronously update all memAgents across all cells
 * 
 * Grow/retract filopodia and lamella veil advance
 * 
 * Activate receptors from local ligand levels
 **/
void World::updateMemAgents(void)
{
    // Init local parameters
    int upto;
    int i, j;

    MemAgent * memp;
    int uptoE = ECagents.size();
    int uptoN, uptoS, uptoSu;
    bool tipDeleteFlag;
    float randomChance;

    bool deleted;

    //--------------------------------------------------------------------------------------------
    // Reset agents and count again
    JunctionAgents.clear(); //used where?
    ALLmemAgents.clear();
    /// Loop over all EC agent and store each mem agent in ALLmemAgents
    for (i = 0; i < uptoE; i++)
    {
        // 
        uptoN = ECagents[i]->nodeAgents.size();
        uptoS = ECagents[i]->springAgents.size();
        uptoSu = ECagents[i]->surfaceAgents.size();
        // Update ALLmemagents
        for (j = 0; j < uptoN; j++) ALLmemAgents.push_back(ECagents[i]->nodeAgents[j]);
        for (j = 0; j < uptoS; j++) ALLmemAgents.push_back(ECagents[i]->springAgents[j]);
        for (j = 0; j < uptoSu; j++) ALLmemAgents.push_back(ECagents[i]->surfaceAgents[j]);

    }
    /// Update the number of mem agents
    upto = ALLmemAgents.size();
    //------------------------------------------------------------------------------------------------------

    //------------------------------------------------------------------------------------------------------
    // Reorder agents randomly
    new_random_shuffle(ALLmemAgents.begin(), ALLmemAgents.end());
    //------------------------------------------------------------------------------------------------------

    //------------------------------------------------------------------------------------------------------
    // Loop over each mem agent: update its prot levels and try to extend/retract filopodia/lamellapodia
    for (i = 0; i < upto; i++)
    {
        // Init temporary variables
        tipDeleteFlag = false;

        memp = ALLmemAgents[i];
        memp->assessed = true;
        memp->addedJunctionList = false;

        // Delete spring agents sitting along filopodia scheduled for deletion during previous fil retraction
        deleted = delete_if_spring_agent_on_a_retracted_fil(memp);

        // Check for deletion -> only for non-deleted mem agents
        if (deleted == false)
        {
            // Reset memAgents active Notch level -> ready for new binding
            memp->activeNotch = 0.0f;
            
            //LC// Maybe keep those options in case code is reused
            //this is needed to tell if triangle positions have changed
            //on the fly surface agent coverage code
            if (on_the_fly_surface_agents == true) memp->store_previous_triangle_pos();

            // Assess local Moore neighbourhood and store data (includes diagonal neighs)
            memp->checkNeighs(false);

            // Determine if agent is on a junction for junctional behaviours
            memp->JunctionTest(true);

            // If the memAgent resides at the tip of a filopodium: check for retraction
            if (memp->FIL == TIP)
            {
                // Generate random prob/chance
                randomChance = new_rand() / (float) NEW_RAND_MAX;

                //LC// Maybe keep those options in case code is reused
                //veil advance for cell migration------------------------
                if (VEIL_ADVANCE == true) {
                    if ((memp->form_filopodia_contact() == true) || (randomChance < RAND_VEIL_ADVANCE_CHANCE)){

                        if ((ANALYSIS_HYSTERESIS == true)&&(memp->Cell != ECagents[0])&&(memp->Cell != ECagents[ECELLS - 1]))
                            memp->veilAdvance();
                        else if(ANALYSIS_HYSTERESIS==false) memp->veilAdvance();
                    }
                }
                //------------------------------------

                // Retract fils if inactive (MAYBE VERIFY THESE CONDITIONS -----------
                if ( ((RAND_FILRETRACT_CHANCE==-1)&&(memp->filTipTimer > FILTIPMAX)) || ((RAND_FILRETRACT_CHANCE>-1)&&(randomChance < RAND_FILRETRACT_CHANCE)) )
                {
                    if (memp->filRetract() == true)
                    {
                        tipDeleteFlag = true;
                        deleteOldGridRef(memp, true);

                        delete memp;
                    }
                        //NEEDED TO CALC CURRENT ACTIN USAEAGE for limit on fil extension //LC//??
                    else
                    {
                        memp->calcRetractDist();
                    }
                }
                //---------------------------------------------------------------------
            }

            // If mem agent not deleted: update receptor activities and check for extension of fil
            if (tipDeleteFlag == false)
            {
                // Reset active dimers to 
                memp->R2R2active = 0.0f;
                memp->R2R3active = 0.0f;
                memp->R3R3active = 0.0f;

                // Apply VEGFRreponse() //it seems that it's not applied in some cases when doing hysteresis simulation
                if ((ANALYSIS_HYSTERESIS == true)&&(memp->Cell != ECagents[0])&&(memp->Cell != ECagents[ECELLS - 1]))
                {
                    if (memp->vonNeu == true) memp->VEGFRresponse();
                }
                else if(ANALYSIS_HYSTERESIS == false)
                {
                    if (memp->vonNeu == true) memp->VEGFRresponse();
                }

                // Apply NotchResponse() in junction agents
                if (memp->junction == true) memp->NotchResponse();

                // Apply TokenTrading()
                ///pass actin to nearest nodes Agent if a surfaceAgent, or further towards tip nodeagent if in a filopodium; lose all if not active
                if ((ANALYSIS_HYSTERESIS == true)&&(memp->Cell != ECagents[0])&&(memp->Cell != ECagents[ECELLS - 1]))
                {
                    //memp->ActinFlow();
                    memp->TokenTrading();
                }
                else if(ANALYSIS_HYSTERESIS == false)
                {
                      //memp->ActinFlow();
                      memp->TokenTrading();
                }
            }

        }
    }

    // the force of new memAgent movements made in functions above are conveyed through the springs following Hookes Law to move all memAgents within the mesh
    if ((ANALYSIS_HYSTERESIS == true)&&(memp->Cell != ECagents[0])&&(memp->Cell != ECagents[ECELLS - 1]))
    {
        calculateSpringAdjustments();
    }
    else if(ANALYSIS_HYSTERESIS == false)
    {
        calculateSpringAdjustments();
    }

}

//------------------------------------------------------------------------------------------
//LC// EC agents update (called at each timestep)
//------------------------------------------------------------------------------------------
/**
 * synchronously update all cell agents after memAgents have updated
 * 
 * Sum new active receptor levels
 * 
 * Calculate GRN new gene expression 
 * 
 * Re-allocate out new levels to memAgents
 **/
void World::updateECagents(void)
{
    // Init local parameters
    int i, j;
    int upto = ECagents.size();
    vector<float> DLL4_timestep;
    vector<float> VEGFR2_timestep;
    vector<float> VEGFR3_timestep;
    vector<float> activeNotch_timestep;
    vector<float> activeR2R2_timestep;
    vector<float> activeR2R3_timestep;
    vector<float> activeR3R3_timestep;

    // Loop over each EC agent: apply protein and actin updates
    for (j = 0; j < upto; j++)
    {
        if (ANALYSIS_COMS == true) ECagents[j]->store_cell_COM(); //to see cell movement, monitor its centre of mass

        // Update actin used after mem agent updates
        ECagents[j]->calcCurrentActinUsed(); //determine overall actin level after filopodia dynamics in memAgent update

        // Update protein values after mem agent updates (+ d)
        ECagents[j]->updateProteinTotals(); //total up the memAgents new active receptor levels, add to time delay stacks 

        // Apply GRN rules on updates protein values
        ECagents[j]->GRN(); //use the time delayed active receptor levels (time to get to get to nucleus+transcription) to calculate gene expression changes

        // Deal with spring integrity etc (keep it like that)
        ECagents[j]->newNodes(); //add new nodes or delete them if springs size is too long/too short (as filopodia have nodes and adhesions along them at 2 micron intervals

        DLL4_timestep.push_back(ECagents[j]->Dll4tot);
        VEGFR2_timestep.push_back(ECagents[j]->VEGFR2tot);
        VEGFR3_timestep.push_back(ECagents[j]->VEGFR3tot);
        activeNotch_timestep.push_back(ECagents[j]->activeNotchtot);
        activeR2R2_timestep.push_back(ECagents[j]->activeR2R2tot);
        activeR2R3_timestep.push_back(ECagents[j]->activeR2R3tot);
        activeR3R3_timestep.push_back(ECagents[j]->activeR3R3tot);
        // DLL4_evolution_file << ECagents[j]->Dll4tot << ", ";
        // VEGFR2_evolution_file << ECagents[j]->VEGFR2tot << ", ";
        // VEGFR3_evolution_file << ECagents[j]->VEGFR3tot << ", ";
        // activeNotch_evolution_file << ECagents[j]->activeNotchtot << ", ";
        // activeR2R2_evolution_file << ECagents[j]->activeR2R2tot << ", ";
        // activeR2R3_evolution_file << ECagents[j]->activeR2R3tot << ", ";
        // activeR3R3_evolution_file << ECagents[j]->activeR3R3tot << ", ";
    }

    if (timeStep < MAXtime)
    {
        std::ostream_iterator<float> output_iterator(RUNSfile, ",");
        std::copy(DLL4_timestep.begin(), DLL4_timestep.end(), output_iterator);
        std::copy(VEGFR2_timestep.begin(), VEGFR2_timestep.end(), output_iterator);
        std::copy(VEGFR3_timestep.begin(), VEGFR3_timestep.end(), output_iterator);
        std::copy(activeNotch_timestep.begin(), activeNotch_timestep.end(), output_iterator);
        std::copy(activeR2R2_timestep.begin(), activeR2R2_timestep.end(), output_iterator);
        std::copy(activeR2R3_timestep.begin(), activeR2R3_timestep.end(), output_iterator);
        std::copy(activeR3R3_timestep.begin(), activeR3R3_timestep.end(), output_iterator);
    }
    //RUNSfile << endl;
    // DLL4_evolution_file << endl;
    // VEGFR2_evolution_file << endl;
    // VEGFR3_evolution_file << endl;
    // activeNotch_evolution_file << endl;
    // activeR2R2_evolution_file << endl;
    // activeR2R3_evolution_file << endl;
    // activeR3R3_evolution_file << endl;

    // Loop over each EC agent: update spring agents
    for (j = 0; j < (int) ECagents.size(); j++)
    {
        // clear all spring agents from previous time step - all springs have been adjusted so need new arrangement of spring agents
        ECagents[j]->removeSpringAgents();

    }

    // Loop over each EC agent: update grid agents (?)
    for (j = 0; j < (int) ECagents.size(); j++)
    {
        //voxellise new spring and surface positions of mesh
        ECagents[j]->gridAgents();

        //LC// Maybe keep those options in case code is reused
        //faster way to do it for debugging versions, but not correct to use in main simulations!!!
        if (on_the_fly_surface_agents == true) ECagents[j]->remove_DoubledUp_SurfaceAgents();
    }

    // Loop over each EC agent: allocate proteins to mem agents
    for (j = 0; j < (int) ECagents.size(); j++) {
        //distribute back out the new VR-2 and Dll4 and Notch levels to voxelised memAgents across the whole new cell surface.
        ECagents[j]->allocateProts();

        //LC// Maybe keep those options in case code is reused
        //use analysis method in JTB paper to obtain tip cell numbers, stability of S&P pattern etc. requird 1 cell per cross section in vessel (PLos/JTB cell setup)
        if (ANALYSIS_JTB_SP_PATTERN == true) ECagents[j]->calcStability();
    }
}


//------------------------------------------------------------------------------------------
//LC// Calculate spring forces and move memAgents in the mesh (called at each timestep in updateMemAgent())
//------------------------------------------------------------------------------------------
/**
 * Calculate spring forces and move memAgents in the mesh
 * 
 * Goes asynchronously currently in the same order through all memAgents in the cells. could randomise in future.
 **/
void World::calculateSpringAdjustments(void)
{
    // Init local variables
    int i, j, upto, uptoE;
    uptoE = ECagents.size();
    MemAgent* memp;

    // Loop over each EC: update mem agents nodes
    for (j = 0; j < uptoE; j++)
    {
        upto = ECagents[j]->nodeAgents.size(); //number of mem agents in cell
        // Loop over each mem agent: update mem agents nodes
        for (i = 0; i < upto; i++)
        {
            memp = ECagents[j]->nodeAgents[i];

            if (memp->FA != true) memp->calcForce();
        }
    }
}


//------------------------------------------------------------------------------------------
//LC// Boom explosion
//------------------------------------------------------------------------------------------
/**
 * delete all agents and data stored if running multiple runs on a cluster (called in destructor) avoids memory leaks :)
 * 
 **/
void World::destroyWorld(void)
{
    // Init local variables
    EC* ec;
    MemAgent* mp;
    Filopodia* fp;
    Spring* sp;
    Contact* cp;
    Macrophage* map;
    int i, j, k;

    // Loop over the whole grid (x, y, z)
    for (i = 0; i < xMAX; i++)
    {
        for (j = 0; j < yMAX; j++)
        {
            for (k = 0; k < zMAX; k++)
            {
                if (grid[i][j][k].type == E)
                {
                    if (grid[i][j][k].Eid != NULL) delete grid[i][j][k].Eid;
                }
            }
        }
    }

    // Loop over EC agents until there are not any left
    while (!ECagents.empty())
    {
        ec = ECagents.back();

        // Loop over EC's mem agents (node) until there are not any left
        while (!ec->nodeAgents.empty()) 
        {
            mp = ec->nodeAgents.back();
            ec->nodeAgents.pop_back();

            delete mp;
        }

        // Loop over EC's mem agents (spring) until there are not any left
        while (!ec->springAgents.empty())
        {
            mp = ec->springAgents.back();
            ec->springAgents.pop_back();
            delete mp;
        }

        // Loop over EC's mem agents (surface) until there are not any left
        while (!ec->surfaceAgents.empty())
        {
            mp = ec->surfaceAgents.back();
            ec->surfaceAgents.pop_back();
            delete mp;
        }

        // Loop over EC's spring (?) until there are not any left
        while (!ec->Springs.empty()) {
            sp = ec->Springs.back();
            ec->Springs.pop_back();
            delete sp;
        }

        // Remove EC agent from list
        ECagents.pop_back();
        delete ec;
    }

    while (!filopodia.empty()) {

        fp = filopodia.back();
        filopodia.pop_back();
        delete fp;
    }

    while (!contacts.empty()) {

        cp = contacts.back();
        contacts.pop_back();
        delete cp;
    }
    while (!macrophages.empty()) {

        map = macrophages.back();
        macrophages.pop_back();
        delete map;
    }

}


//------------------------------------------------------------------------------------------
//LC// Make movie from screenshots
//------------------------------------------------------------------------------------------
/**
 * 
 * @param movie
 * 
 * Save tiff frames and make movie using imageMagick
 * 
 * !function needs updating! currently using imageMagick, but use something better
 **/
void World::movieMaking(int movie)
{
    //movie making
    if (movie == 1) system("rm -vf movie/*");
    if (movie == 2)
    {
        sprintf(fname, "movie/frame%05i.tga", timeStep);
        fprintf(stdout, "%s %i\n", fname, timeStep);

        //---------------------------
        //uncomment this function (defined in display.cpp) if you can get imageMAgick to run properly then this will actually save the tiff frames in order to make the movie from them
        //SaveFrame(fname);
        //---------------------------

    }
    if (movie == 3)
    {
        //576 x 620 - width must be divisible by 8 when using this imageMagick method!!
        system("cd movie ; /opt/local/bin/mencoder mf:// -mf w=800:h=800:fps=25:type=tga -ovc lavc -lavcopts vcodec=mpeg1video:keyint=25:vbitrate=3000:trell:mbd=2  -oac copy -o output.mpg");
        system(fname);


    }
}


//------------------------------------------------------------------------------------------
//LC// Scale protein levels to cell size (called once, at World creation)
//------------------------------------------------------------------------------------------
/**
 * Scale proteins for coarser grid - currently only 0.25 scaling of prots used (as surface area of cell reduced by this when grid changed from half micron cubes (3D grid site) to one micron cubes
 **/
void World::scale_ProtLevels_to_CellSize(void)
{
    // Init local variables
    float Scale;

    if (CELL_SETUP == 1)
    {
        Scale = 100.0f;
        actinMax = 512;
    } 
    else if ((ECcross == 2)&&(ECwidth == 20) && (vesselRadius == 6)) Scale = 48.6f;
    else
    {
        Scale = 25;
        actinMax = 128;
    }

    // DEFINE INITIAL NUMBER OF PROTEINS (+ limits)
    //half the number of receptors for two cells per vessel cross section - actually 48% as this is the exact ratio of memAgents when two cell.
    if (ECcross == 2)
    {
        NotchNorm = (10000.0f / 100.0f) * Scale;

        MAX_dll4 = (10000.0f / 100.0f) * Scale;
        
        // if (lymphatic==0)
        // {
        //     VEGFR2NORM = (2 * 31714.0f / 100.0f) * Scale; //total of receptors it will maintain if all else is equal - divides out to M agents
        //     // VEGFR3NORM = VEGFR2NORM * 0.58;
        //     VEGFR3NORM = VEGFR2NORM;
        //     cout << "bien vu";
        // }
        // else if (lymphatic==1)
        // {
        //     VEGFR2NORM = 1.3 * (2 * 31714.0f / 100.0f) * Scale;
        //     VEGFR3NORM = 1.153846 * VEGFR2NORM;
        //     cout << "mal vu";
        // }
        
        //LC// multiply by 2 for dimerization so that R2R2 gives the same result as VEGFR2 alone in previous version of code
        VEGFR2NORM = (2 * 31714.0f / 100.0f) * Scale; //total of receptors it will maintain if all else is equal - divides out to M agents
        VEGFR3NORM = VEGFR2NORM * VEGFR3scale;
        //LC// multiply by 2 for dimerization so that R2R2 gives the same result as VEGFR2 alone in previous version of code
        VEGFR2min = (2 * 689.0f / 100.0f) * Scale;
        VEGFR3min = (2 * 689.0f / 100.0f) * Scale;
    }
    else
    {
        NotchNorm = 10000.0f;

        MAX_dll4 = 10000.0f;

        if (lymphatic==0)
        {
            VEGFR2NORM = (2 * 31714.0f); //total of receptors it will maintain if all else is equal - divides out to M agents
            VEGFR3NORM = VEGFR2NORM * 0.58;
            // cout << "JE RENTRE DANS LE IF BORDEL - 0" << endl << endl << endl;
        }
        else if (lymphatic==1)
        {
            VEGFR2NORM = 1.3 * (2 * 31714.0f);
            VEGFR3NORM = 1.153846 * VEGFR2NORM;
            // cout << "JE RENTRE DANS LE IF BORDEL - 1" << endl << endl << endl;
        }
        else if (lymphatic==2)
        {
            VEGFR2NORM = 2 * 31714.0f;
            VEGFR3NORM = VEGFR2NORM * VEGFR3scale;
        }

        // //LC// multiply by 2 for dimerization so that R2R2 gives the same result as VEGFR2 alone in previous version of code
        // VEGFR2NORM = 2 * 31714.0f; //scaled to fit with first model - so each memagent has same number of recs - new arrangment means diff number of initial memagents
        // VEGFR3NORM = VEGFR2NORM * VEGFR3scale;
        //LC// multiply by 2 for dimerization so that R2R2 gives the same result as VEGFR2 alone in previous version of code
        VEGFR2min = 2 * 689.0f;
        VEGFR3min = 2 * 689.0f;
    }
}


//------------------------------------------------------------------------------------------
//LC// Tell if spring agent will delete or not (called at every timestep, for each mem agent)
//------------------------------------------------------------------------------------------
/**
 * 
 * @param memp
 * @return bool flag deleted
 * 
 * in filretract() spring agents are scheduled and flagged for deletion if the filopodia spring they belong to has been retracted back
 * 
 * this function actually deletes them from memory
 */
bool World::delete_if_spring_agent_on_a_retracted_fil(MemAgent* memp) 
{
    // Init local variables
    int k = 0;
    int flag = 0;
    bool deleted = false;
    vector<MemAgent*>::iterator L;

    if ((memp->node == false) && (memp->deleteFlag == true))
    {
        deleted = true;
        //remove the tip node from cells list
        k = 0;
        flag = 0;
        do
        {
            if (memp->Cell->springAgents[k] == memp)
            {
                flag = 1;
                L = memp->Cell->springAgents.begin();
                memp->Cell->springAgents.erase(L + k);
            }
            k++;
        } while ((k < (int) memp->Cell->springAgents.size()) && (flag == 0));

        if (flag == 0)
        {
            cout << "deleting spring agent in main: hasnt found in list" << endl;
            cout.flush();
        }
        delete memp;
    }
    return (deleted);
}


//------------------------------------------------------------------------------------------
//LC// ?? (seems to be called only when GRAPHICS==true)
//------------------------------------------------------------------------------------------
void World::store_normals(void)
{
    Coordinates cross;

    //top face z+L for all do minus sign of this for z face

    cross.x = 0.0;
    cross.y = 0.0;
    cross.z = 1.0;
    store_cube_normals.push_back(cross);

    cross.x = 1.0;
    cross.y = 0.0;
    cross.z = 0.0;
    //left face x+L for all

    store_cube_normals.push_back(cross);

    cross.x = 0.0;
    cross.y = 1.0;
    cross.z = 0.0;
    //back face y+L for all

    store_cube_normals.push_back(cross);
}


//------------------------------------------------------------------------------------------
//LC// Random number generator
//------------------------------------------------------------------------------------------
int new_rand()
{
    return (int)dist(g);
}

//------------------------------------------------------------------------------------------
//LC// Random number generator
//------------------------------------------------------------------------------------------
float new_rand_beta(float alpha, float beta)
{
    float X, Y;
    std::gamma_distribution<double> dist_X(alpha, 1.0);
    std::gamma_distribution<double> dist_Y(beta, 1.0);
    X = dist_X(g);
    Y = dist_Y(g);
    return (float)(X/(X+Y));
}

//------------------------------------------------------------------------------------------
//LC// Create statistics file
//------------------------------------------------------------------------------------------
void create_statistics_file(string statisticsFilename)
{
	ofstream statisticsFile(statisticsFilename);
	statisticsFile.close();
}


//------------------------------------------------------------------------------------------
//LC// Write in statistics file
//------------------------------------------------------------------------------------------
void write_to_statistics_file(string statisticsFilename, string line)
{
	ofstream statisticsFile;
	statisticsFile.open(statisticsFilename, std::ios_base::app);
	if (statisticsFile.is_open()) statisticsFile << line;
	statisticsFile.close();
}


//------------------------------------------------------------------------------------------
//LC// Get current time when starting and finishing simulation
//------------------------------------------------------------------------------------------
std::time_t get_current_time() {
	auto time = std::chrono::system_clock::now();
	std::time_t current_time = std::chrono::system_clock::to_time_t(time);
	return current_time;
}


//------------------------------------------------------------------------------------------
//LC// Associate string to time (start or end)
//------------------------------------------------------------------------------------------
std::string format_time_string(std::time_t time, bool start) {
	// N.B. Function should be called at start and end of simulation only, for logging purposes.
	string time_entry, time_string;

	if (start)
    {
		time_entry = "Start time,";
	}
    else
    {
		time_entry = "End time,";
	}

	time_string = std::ctime(&time);
	time_entry = time_entry + time_string;

	return time_entry;
}
