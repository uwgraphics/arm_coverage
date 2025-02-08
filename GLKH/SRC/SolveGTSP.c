#include "LKH.h"
#include <unistd.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/inotify.h>
#include <fcntl.h>
#include <errno.h>

GainType SolveTSP(int Dimension, char *ParFileName,
                  char *TourFileName, int *Tour, GainType Optimum,
                  GainType Deduction);
static GainType mySolveTSP(int Dimension, char *ParFileName,
                  char *TourFileName, int *Tour, GainType Optimum,
                  GainType Deduction, int Clusters);

enum TourType { INITIAL, INPUT, MERGE, SUBPROBLEM };
static void WriteFullTour(enum TourType Type, int Dimension,
                          char *TourFileName, int Id);
static void ExtractTour(int *Tour, int *GTour, int Clusters);
static void ReadTourFile(char *TourFileName, int *Tour, int Dimension);

#define EVENT_SIZE  ( sizeof (struct inotify_event) )
#define BUF_LEN     ( 1024 * ( EVENT_SIZE + 16 ) )

/*
 * The SolveGTSP function solves an E-GTSP instance. The resulting g-tour is
 * returned in the parameter GTour, and its cost as the value of the 
 * function.
 *
 * The algorithm is as follows:
 *. 1. Transform the E-GTSP instance into an asymmetric TSP instance.
 *  2. Write the TSP instance to a problem file.
 *  3. Write suitable parameter values to a parameter file.
 *  4. Execute LKH given these two files (by calling solveTSP).
 *  5. Extract the g-tour from the TSP solution tour by picking the 
 *     first vertex from each cluster in the TSP tour.
 */

GainType SolveGTSP(int *GTour)
{
    int i, j, Dist, Clusters = 0;
    Cluster *Cl;
    Node *From, *To;
    FILE *ParFile, *ProblemFile;
    char ParFileName[256], ProblemFileName[256], TourFileName[256],
        NewInitialTourFileName[256] = { 0 }, 
        NewInputTourFileName[256] = { 0 }, 
        NewSubproblemTourFileName[256] = { 0 },
        **NewMergeTourFileName, 
        Prefix[256];
    GainType M, Cost;
    int *Tour, *Used;

    assert(NewMergeTourFileName =
           (char **) malloc(MergeTourFiles * sizeof(char *)));
    for (i = 0; i < MergeTourFiles; i++)
        assert(NewMergeTourFileName[i] = (char *) malloc(256));
    for (Cl = FirstCluster; Cl; Cl = Cl->Next) {
        Clusters++;
        From = Cl->First;
        do
            From->V = Clusters;
        while ((From = From->Next) != Cl->First);
    }
    assert(Clusters == GTSPSets);

    M = Clusters < 2 ? 0 : INT_MAX / 4 / Precision;

    sprintf(Prefix, "%s.pid%d", Name, getpid());

    /* Create the problem file */
    sprintf(ProblemFileName, "TMP/%s.atsp", Prefix);
    assert(ProblemFile = fopen(ProblemFileName, "w"));
    fprintf(ProblemFile, "NAME : %s.gtsp\n", Prefix);
    fprintf(ProblemFile, "TYPE : ATSP\n");
    if (ProblemType != ATSP)
        fprintf(ProblemFile, "DIMENSION : %d\n", Dimension);
    else
        fprintf(ProblemFile, "DIMENSION : %d\n", DimensionSaved);

    if (1) {
        fprintf(ProblemFile, "EDGE_DATA_FORMAT : EDGE_LIST\n");
        fprintf(ProblemFile, "EDGE_DATA_SECTION: \n");
        /* Transform the GTSP into an ATSP */
        for (i = 1; i <= DimensionSaved; i++) {
            From = &NodeSet[i];
            for (j = 1; j <= DimensionSaved; j++) {
                if (i == j)
                    continue;
                else {
                    To = &NodeSet[j];
                    Dist = To == From->Next ? 0 :
                        To->V == From->V ? M * 2:
                        Distance(From->Next, To) + M;
                        // Distance(From->Next, To) < 536870911 ? Distance(From->Next, To) + M : -1;
                        // (ProblemType != ATSP ? Distance(From->Next, To) :
                        //  From->Next->C[j]);
                    if (Dist != -1 && (Dist * Precision < 0 ||
                        Dist * Precision / Precision != Dist))
                        printff("*** PRECISION (= %d) is too large. Dist: %d", Precision, Dist);
                    // fprintf(ProblemFile, "%d ", Dist);
                    if (Dist < M * 2)
                        fprintf(ProblemFile, "%d %d %d\n", i, j, Dist);
                }
            }
        }
    } else {
        fprintf(ProblemFile, "EDGE_WEIGHT_TYPE : EXPLICIT\n");
        fprintf(ProblemFile, "EDGE_WEIGHT_FORMAT : FULL_MATRIX\n");
        fprintf(ProblemFile, "EDGE_WEIGHT_SECTION\n");
        for (i = 1; i <= DimensionSaved; i++) {
            From = &NodeSet[i];
            for (j = 1; j <= DimensionSaved; j++) {
                if (i == j)
                    fprintf(ProblemFile, "999999 ");
                else {
                    To = &NodeSet[j];
                    Dist = To == From->Next ? 0 :
                        To->V == From->V ? M * 2:
                        Distance(From->Next, To) + M;
                        // (ProblemType != ATSP ? Distance(From->Next, To) :
                        //  From->Next->C[j]);
                    if (Dist * Precision < 0 ||
                        Dist * Precision / Precision != Dist)
                        printff("*** PRECISION (= %d) is too large", Precision);
                    fprintf(ProblemFile, "%d ", Dist);
                }
            }
            fprintf(ProblemFile, "\n");
        }
    }

    fprintf(ProblemFile, "EOF\n");
    fclose(ProblemFile);

    /* Create the parameter file */
    sprintf(ParFileName, "TMP/%s.par", Prefix);
    assert(ParFile = fopen(ParFileName, "w"));
    fprintf(ParFile, "PROBLEM_FILE = TMP/%s.atsp\n", Prefix);
    fprintf(ParFile, "ASCENT_CANDIDATES = %d\n", AscentCandidates);
    fprintf(ParFile, "BACKBONE_TRIALS = %d\n", BackboneTrials);
    if (Backtracking)
        fprintf(ParFile, "BACKTRACKING  = YES\n");
    for (i = 0; i < CandidateFiles; i++)
        fprintf(ParFile, "CANDIDATE_FILE = %s\n", CandidateFileName[i]);
    fprintf(ParFile, "CANDIDATE_SET_TYPE = %s\n",
            CandidateSetType == ALPHA ? "ALPHA" :
            CandidateSetType == POPMUSIC ? "POPMUSIC" : "");
    if (Excess > 0)
        fprintf(ParFile, "EXCESS = %g\n", Excess);
    if (!Gain23Used)
        fprintf(ParFile, "GAIN23 = NO\n");
    if (!GainCriterionUsed)
        fprintf(ParFile, "GAIN_CRITERION = NO\n");
    fprintf(ParFile, "INITIAL_PERIOD = %d\n", InitialPeriod);
    if (InitialTourAlgorithm != WALK)
        fprintf(ParFile, "INITIAL_TOUR_ALGORITHM = %s\n",
                InitialTourAlgorithm ==
                NEAREST_NEIGHBOR ? "NEAREST-NEIGHBOR" :
                InitialTourAlgorithm == GREEDY ? "GREEDY" : "");
    fprintf(ParFile, "INITIAL_STEP_SIZE = %d\n", InitialStepSize);
    if (InitialTourFileName) {
        sprintf(NewInitialTourFileName, "TMP/%s.initial.tour", Prefix);
        WriteFullTour(INITIAL, DimensionSaved, NewInitialTourFileName, 0);
        fprintf(ParFile, "INITIAL_TOUR_FILE = %s\n",
                NewInitialTourFileName);
    }
    fprintf(ParFile, "INITIAL_TOUR_FRACTION = %0.3f\n",
            InitialTourFraction);
    if (InputTourFileName) {
        sprintf(NewInputTourFileName, "TMP/%s.input.tour", Prefix);
        WriteFullTour(INPUT, DimensionSaved, NewInputTourFileName, 0);
        fprintf(ParFile, "INPUT_TOUR_FILE = %s\n", NewInputTourFileName);
    }
    fprintf(ParFile, "KICK_TYPE = %d\n", KickType);
    fprintf(ParFile, "MAX_BREADTH = %d\n", MaxBreadth);
    fprintf(ParFile, "MAX_CANDIDATES = %d%s\n", MaxCandidates,
            CandidateSetSymmetric ? " SYMMETRIC" : "");
    fprintf(ParFile, "MAX_SWAPS = %d\n", MaxSwaps);
    fprintf(ParFile, "MAX_TRIALS = %d\n", MaxTrials);
    for (i = 0; i < MergeTourFiles; i++) {
        sprintf(NewMergeTourFileName[i],
                "TMP/%s.merge%d.tour", Prefix, i + 1);
        WriteFullTour(MERGE, DimensionSaved, NewMergeTourFileName[i], i);
        fprintf(ParFile, "MERGE_TOUR_FILE = %s\n",
                NewMergeTourFileName[i]);
    }
    fprintf(ParFile, "MOVE_TYPE = %d\n", MoveType);
    if (NonsequentialMoveType >= 4)
        fprintf(ParFile, "NONSEQUENTIAL_MOVE_TYPE = %d\n",
                NonsequentialMoveType);
    if (Optimum != MINUS_INFINITY)
        fprintf(ParFile, "OPTIMUM = " GainFormat "\n",
                Optimum + Clusters * M);
    fprintf(ParFile, "PATCHING_A = %d %s\n", PatchingA,
            PatchingARestricted ? "RESTRICTED" :
            PatchingAExtended ? "EXTENDED" : "");
    fprintf(ParFile, "PATCHING_C = %d %s\n", PatchingC,
            PatchingCRestricted ? "RESTRICTED" :
            PatchingCExtended ? "EXTENDED" : "");
    if (PiFileName)
        fprintf(ParFile, "PI_FILE = %s\n", PiFileName);
    fprintf(ParFile, "POPMUSIC_INITIAL_TOUR = %s\n",
            POPMUSIC_InitialTour ? "YES" : "NO");
    fprintf(ParFile, "POPMUSIC_MAX_NEIGHBORS = %d\n", POPMUSIC_MaxNeighbors);
    fprintf(ParFile, "POPMUSIC_SAMPLE_SIZE = %d\n", POPMUSIC_SampleSize);
    fprintf(ParFile, "POPMUSIC_SOLUTIONS = %d\n", POPMUSIC_Solutions);
    fprintf(ParFile, "POPMUSIC_TRIALS = %d\n", POPMUSIC_Trials);
    fprintf(ParFile, "POPULATION_SIZE = %d\n", MaxPopulationSize);
    fprintf(ParFile, "PRECISION = %d\n", Precision);
    fprintf(ParFile, "RECOMBINATION = %s\n",
                     Recombination == GPX2 ? "GPX2" : "IPT");
    if (!RestrictedSearch)
        fprintf(ParFile, "RESTRICTED_SEARCH = NO\n");
    fprintf(ParFile, "RUNS = %d\n", Runs);
    fprintf(ParFile, "SEED = %d\n", Seed);
    if (!StopAtOptimum)
        fprintf(ParFile, "STOP_AT_OPTIMUM = NO\n");
    if (!Subgradient)
        fprintf(ParFile, "SUBGRADIENT = NO\n");
    if (SubproblemSize > 0)
        fprintf(ParFile, "SUBPROBLEM_SIZE = %d\n", 2 * SubproblemSize);
    if (SubproblemTourFileName) {
        sprintf(NewSubproblemTourFileName, "TMP/%s.subproblem.tour",
                Prefix);
        WriteFullTour(SUBPROBLEM, DimensionSaved, NewSubproblemTourFileName, 0);
        fprintf(ParFile, "SUBPROBLEM_TOUR_FILE = %s\n", 
                NewSubproblemTourFileName);
    }
    fprintf(ParFile, "SUBSEQUENT_MOVE_TYPE = %d\n", SubsequentMoveType);
    if (!SubsequentPatching)
        fprintf(ParFile, "SUBSEQUENT_PATCHING = NO\n");
    if (TimeLimit != DBL_MAX)
        fprintf(ParFile, "TIME_LIMIT = %0.1f\n", TimeLimit);
    sprintf(TourFileName, "TMP/%s.tour", Prefix);
    fprintf(ParFile, "OUTPUT_TOUR_FILE = %s\n", TourFileName);
    fprintf(ParFile, "TRACE_LEVEL = %d\n",
            TraceLevel == 0 ? 1 : TraceLevel);
    fclose(ParFile);

    /* Solve the ATSP */
    assert(Tour = (int *) malloc((DimensionSaved + 1) * sizeof(int)));
    Cost =
        mySolveTSP(DimensionSaved, ParFileName, TourFileName,
                 Tour, Optimum, Clusters * M, Clusters);
    // unlink(ParFileName);
    // unlink(ProblemFileName);
    unlink(NewInitialTourFileName);
    unlink(NewInputTourFileName);
    for (i = 0; i < MergeTourFiles; i++)
        unlink(NewMergeTourFileName[i]);
    unlink(NewSubproblemTourFileName);

    // ExtractTour(Tour, GTour, Clusters);
    free(Tour);
    return Cost;
}

static GainType mySolveTSP(int Dimension, char *ParFileName,
                  char *TourFileName, int *Tour, GainType Optimum,
                  GainType Deduction, int Clusters)
{
    GainType Cost;

    int *GTour;
    assert(GTour = (int *) malloc((GTSPSets + 1) * sizeof(int)));

    FILE *p, *TourFile;
    int i;
    char Command[256], Key[256], Buffer[256], *Line, *Keyword;
    char Delimiters[] = " :=\n\t\r\f\v\xef\xbb\xbf";

    sprintf(Command, "./LKH %s", ParFileName);
    assert(p = popen(Command, "r"));
    Cost = PLUS_INFINITY;


    int length, j = 0;
    int fd;
    int wd;
    char buffer[BUF_LEN];
    fd = inotify_init();

    if (fd < 0) {
        perror("inotify_init");
        exit(EXIT_FAILURE);
    }

    int flags = fcntl(fd, F_GETFL, 0);
    fcntl(fd, F_SETFL, flags | O_NONBLOCK);

    creat(TourFileName, 0666);

    wd = inotify_add_watch(fd, TourFileName, IN_CLOSE_WRITE);

    // Check if watch was added successfully
    if (wd == -1) {
        perror("inotify_add_watch");
        close(fd);
        exit(EXIT_FAILURE);
    }

    // length = read(fd, buffer, sizeof(buffer));

    // printff("The file %s was modified.\n", TourFileName);

    // if ( length < 0 ) {
    //     perror( "read" );
    // }

    // while ( j < length ) {
    //     struct inotify_event *event = ( struct inotify_event * ) &buffer[ j ];
    //     if (event->len) {
    //         if (event->mask & IN_CREATE) {
    //             printf("The file %s was created.\n", event->name);
    //         } else if (event->mask & IN_DELETE) {
    //             printf("The file %s was deleted.\n", event->name);
    //         } else if (event->mask & IN_MODIFY) {
    //             printf("The file %s was modified.\n", event->name);
    //         }
    //     }
    //     j += EVENT_SIZE + event->len;
    // }

    // ( void ) inotify_rm_watch( fd, wd );
    // ( void ) close( fd );

    GainType LocalCost;
    while (fgets(Buffer, sizeof(Buffer), p)) {
        if (TraceLevel > 1)
            printff("%s", Buffer);
        if (sscanf(Buffer, "%s", Key) > 0) {
            int LocalTrial;
            double LocalTime;
            if (!strcmp(Key, "Cost.min")) {
                char *cp = strchr(Buffer + strlen(Key), '=');
                sscanf(cp + 1, GainInputFormat, &Cost);
                Cost -= Deduction;
            } else if (TraceLevel > 0) { 
                if (!strcmp(Key, "Run")) {
                    char *cp = Buffer + strlen(Key);
                    sscanf(cp + 1, "%d", &Run);
                    sscanf(cp + 1, "%d", &LocalTrial);
                    cp = strchr(cp + 1, '=');
                    sscanf(cp + 1, GainInputFormat, &LocalCost);
                    LocalCost -= Deduction;
                    cp = strchr(cp + 1, 'T');
                    cp = strchr(cp + 1, '=');
                    sscanf(cp + 1, "%lf", &LocalTime);
                    printff("Run %d: Cost = " GainFormat ", ",
                            Run, LocalCost);
                    if (Optimum != MINUS_INFINITY && Optimum != 0)
                        printf("Gap = %0.4f%%, ",
                               100.0 * (LocalCost - Optimum) / Optimum);
                    printff("Time = %0.2f sec.\n\n", LocalTime);
                } else if (!strcmp(Key, "*")) {
                    char *cp = Buffer + strlen(Key);
                    if (sscanf(cp + 1, "%d", &LocalTrial) > 0) {
                        cp = strchr(cp + 1, '=');
                        sscanf(cp + 1, GainInputFormat, &LocalCost);
                        LocalCost -= Deduction;
                        cp = strchr(cp + 1, 'T');
                        cp = strchr(cp + 1, '=');
                        sscanf(cp + 1, "%lf", &LocalTime);
                        printff("# %d: Cost = " GainFormat ", ",
                                LocalTrial, LocalCost);
                        if (Optimum != MINUS_INFINITY && Optimum != 0)
                            printf("Gap = %0.4f%%, ",
                                   100.0 * (LocalCost - Optimum) / Optimum);
                        printff("Time = %0.2f sec.\n", LocalTime);
                    }
                }
            } else if (!strcmp(Key, "Run")) {
                char *cp = Buffer + strlen(Key);
                sscanf(cp + 1, "%d", &Run);
            }

        }

        length = read(fd, buffer, BUF_LEN);

        // Check if read was successful
        if (length < 0) {
            if (errno == EAGAIN) {
                // No data available
                continue;
            } else {
                perror("read");
                close(wd);
                close(fd);
                exit(EXIT_FAILURE);
            }
        }

        // Process each event
        i = 0;
        while (i < length) {
            struct inotify_event *event = (struct inotify_event *) &buffer[i];
            if (!event->len) {
                printff("File %s was modified.\n", TourFileName);
                ReadTourFile(TourFileName, Tour, Dimension);
                // print out tour
                // printff("in loop Tour: ");
                // for (i = 0; i <= Dimension; i++) {
                //     printff("%d ", Tour[i]);
                // }
                // printff("\n");
                ExtractTour(Tour, GTour, Clusters);
                // printout GTour
                printff("GTour: ");
                for (i = 0; i <= Clusters; i++) {
                    printff("%d ", GTour[i]);
                }
                printff("\n");

                printff("LocalCost = " GainFormat "\n", LocalCost);
                
                LocalCost = PostOptimize(GTour, LocalCost);

                printff("PostOptimize GTour: ");
                for (i = 0; i <= Clusters; i++) {
                    printff("%d ", GTour[i]);
                }
                printff("\n");

                WriteTour(OutputTourFileName, GTour, LocalCost);
            }
            i += EVENT_SIZE + event->len;
        }
    }

    pclose(p);
    // ReadTourFile(TourFileName, Tour, Dimension);
    // TourFile = fopen(TourFileName, "r");
    // if (!TourFile)
    //     return PLUS_INFINITY;
    // while ((Line = ReadLine(TourFile))) {
    //     if (!(Keyword = strtok(Line, Delimiters)))
    //         continue;
    //     for (i = 0; i < strlen(Keyword); i++)
    //         Keyword[i] = (char) toupper(Keyword[i]);
    //     if (!strcmp(Keyword, "TOUR_SECTION"))
    //         break;
    // }
    // for (i = 1; i <= Dimension; i++)
    //     fscanf(TourFile, "%d", &Tour[i]);
    // Tour[0] = Tour[Dimension];
    // fclose(TourFile);

    // unlink(TourFileName);
    return Cost;
}

static void ReadTourFile(char *TourFileName, int *Tour, int Dimension)
{
    FILE *TourFile;
    char *Line, *Keyword;
    char Delimiters[] = " :=\n\t\r\f\v\xef\xbb\xbf";
    int i;

    TourFile = fopen(TourFileName, "r");
    if (!TourFile)
        return PLUS_INFINITY;
    while ((Line = ReadLine(TourFile))) {
        if (!(Keyword = strtok(Line, Delimiters)))
            continue;
        for (i = 0; i < strlen(Keyword); i++)
            Keyword[i] = (char) toupper(Keyword[i]);
        if (!strcmp(Keyword, "TOUR_SECTION"))
            break;
    }
    for (i = 1; i <= Dimension; i++)
        fscanf(TourFile, "%d", &Tour[i]);
    Tour[0] = Tour[Dimension];

    // print out tour for debugging
    printff("Read TSP Tour: ");
    for (i = 0; i <= Dimension; i++) {
        printff("%d ", Tour[i]);
    }
    printff("\n");

    fclose(TourFile);
}

static void ExtractTour(int *Tour, int*GTour, int Clusters){
    int i, j;
    int  *Used;
    Node *From, *To;
    /* Extract the g-tour and check it */
    for (i = 1; i <= DimensionSaved; i++) {
        NodeSet[Tour[i - 1]].Suc = &NodeSet[Tour[i]];
    }
    // free(Tour);
    From = FirstNode;
    i = From->V;
    do
        FirstNode = From = From->Suc;
    while (From->V == i);
    assert(Used = (int *) calloc(Clusters + 1, sizeof(int)));
    i = 0;
    do {
        GTour[++i] = From->Id;
        j = From->V;
        if (Used[j])
            eprintf("Illegal g-tour: cluster entered more than once");
        Used[j] = 1;
        while (From->V == j)
            From = From->Suc;
    } while (From != FirstNode);
    free(Used);
    if (i != Clusters)
        eprintf("Illegal g-tour: unvisited cluster(s)");
    GTour[0] = GTour[Clusters];
}

static void WriteFullTour(enum TourType Type, int Dimension,
                         char *TourFileName, int Id)
{
    int *Tour, i = 0;
    Node *N, *From, *To;

    assert(Tour = (int *) malloc((Dimension + 1) * sizeof(int)));
    From = FirstNode;
    while (Type == INITIAL ? !From->InitialSuc :
           Type == INPUT ? !From->InputSuc :
           Type == MERGE ? !From->MergeSuc[Id] :
           Type == SUBPROBLEM ? !From->SubproblemSuc : 0)
        From = From->Suc;
    N = From;
    do {
        To = N;
        do
            Tour[++i] = N->Id;
        while ((N = N->Next) != To);
        N = Type == INITIAL ? N->InitialSuc :
            Type == INPUT ? N->InputSuc :
            Type == MERGE ? N->MergeSuc[Id] :
            Type == SUBPROBLEM ? N->SubproblemSuc : From;
    } while (N != From);
    assert(i == Dimension);
    Tour[0] = Tour[Dimension];
    WriteTour(TourFileName, Tour, -1);
    free(Tour);
}
