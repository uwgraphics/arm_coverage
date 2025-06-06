#include "LKH.h"

/*
 * The WriteTour function writes a tour to file. The tour 
 * is written in TSPLIB format to file FileName. 
 * 
 * The tour is written in "normal form": starting at node 1,
 * and continuing in direction of its lowest numbered
 * neighbor.
 * 
 * Nothing happens if FileName is 0. 
 */

static char *FullName(char *Name, GainType Cost);

void WriteTour(char *FileName, int *Tour, GainType Cost)
{
    FILE *TourFile;
    int i, j, n, Forwards;
    char *FullFileName;
    time_t Now;

    n =  ProblemType != ATSP ? GTSPSets : 2 * GTSPSets;

    printff("WriteTour\n");
    for (i = 1; i <= n; i++)
        printff("%d ", Tour[i]);
    printff("\n");

    if (FileName == 0)
        return;
    FullFileName = FullName(FileName, Cost);
    Now = time(&Now);
    if (TraceLevel >= 1)
        printff("Writing%s: \"%s\" ... ",
                FileName == TourFileName ? " TOUR_FILE" :
                FileName == OutputTourFileName ? " OUTPUT_TOUR_FILE" : "",
                FullFileName);
    assert(TourFile = fopen(FullFileName, "w"));
    fprintf(TourFile, "NAME : %s." GainFormat ".tour\n", Name, Cost);
    fprintf(TourFile, "COMMENT : Length = " GainFormat "\n", Cost);
    fprintf(TourFile, "COMMENT : Found by GLKH [Keld Helsgaun] %s",
            ctime(&Now));
    fprintf(TourFile, "TYPE : TOUR\n");

    fprintf(TourFile, "DIMENSION : %d\n", n);
    fprintf(TourFile, "TOUR_SECTION\n");

    i = 1;
    for (j = 1; j <= n; j++)
        if (Tour[j] < Tour[i])
            i = j;
    Forwards = 1;
    for (j = 1; j <= n; j++) {
        fprintf(TourFile, "%d\n", Tour[i]);
        if (Forwards) {
            if (++i > n)
                i = 1;
        } else if (--i < 1)
            i = n;
    }
    fprintf(TourFile, "-1\nEOF\n");
    fclose(TourFile);
    if (TraceLevel >= 1)
        printff("done\n");
    free(FullFileName);
}

/*
 * The FullName function returns a copy of the string Name where all 
 * occurrences of the character '$' have been replaced by Cost.        
 */

static char *FullName(char *Name, GainType Cost)
{
    char *NewName = 0, *CostBuffer, *Pos;

    if (!(Pos = strstr(Name, "$"))) {
        assert(NewName = (char *) calloc(strlen(Name) + 1, 1));
        strcpy(NewName, Name);
        return NewName;
    }
    assert(CostBuffer = (char *) malloc(400));
    sprintf(CostBuffer, GainFormat, Cost);
    do {
        free(NewName);
        assert(NewName =
               (char *) calloc(strlen(Name) + strlen(CostBuffer) + 1, 1));
        strncpy(NewName, Name, Pos - Name);
        strcat(NewName, CostBuffer);
        strcat(NewName, Pos + 1);
        Name = NewName;
    }
    while ((Pos = strstr(Name, "$")));
    free(CostBuffer);
    return NewName;
}
