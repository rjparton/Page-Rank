// Written by Robert Parton
// 17 November 2022

#include <assert.h>
#include <ctype.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Graph.h"
#include "List.h"

#define MAX_STRLEN 100

List readCollectionFile();
Graph createGraph(List l);
static void insertEdges(Graph g, List l, FILE *fp, char *url, Node curr);
static int getUrlIndex(List l, char *url);
List calculatePageRank(List l, Graph g, double d, double diffPR, 
                        int maxIterations);
static double getPageWeight(List l, Graph g, Graph gWin, Graph gWout, Node pi);
static void initialiseRankAndDegree(List l, Graph g, double N);
static double calculateDiff(List l);
static Graph setGraphWin(List l, Graph g);
static Graph setGraphWout(List l, Graph g);
static double calculateWin(List l, Graph g, Node pj, Node pi);
static double incomingDegree(Graph g, int urlA);
static double calculateWout(List l, Graph g, Node pj, Node pi);
static double outgoingDegree(Graph g, int urlA);
static bool isAdjacent(Graph g, int v, int w);
static Node getNode(List l, int index);

int main(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s dampingFactor diffPR maxIterations\n",
                argv[0]);
        return EXIT_FAILURE;
    }

    // Convert inputs from strings to integers
    double d = atof(argv[1]);
    double diffPR = atof(argv[2]);
    int maxIterations = atoi(argv[3]);

    // Read URLs and store in a Linked List
    List urlList = readCollectionFile();

    // Create Graph Matrix for urlList
    Graph urlGraph = createGraph(urlList);

    // Calculate the page ranks for each url
    urlList = calculatePageRank(urlList, urlGraph, d, diffPR, maxIterations);
    ListSort(urlList);

    ListPrint(urlList);

    ListFree(urlList);
    GraphFree(urlGraph);
    return 0;
}

//
// Helper Functions
//

// Reads the collection.txt file and creates a Linked List containing the URLs
List readCollectionFile() {
    List urlList = ListNew();
    FILE *fp = fopen("collection.txt", "r");
    if (fp == NULL) {
        fprintf(stderr, "fopen\n");
        exit(EXIT_FAILURE);
    }
    
    // Create the linked list to store the urls
    char *url = malloc(MAX_STRLEN * sizeof(char) + 1);
    if (url == NULL) {
        fprintf(stderr, "error: out of memory\n");
        exit(EXIT_FAILURE);
    }

    while (fscanf(fp, "%s", url) != EOF) {
        ListAppend(urlList, url);
    }

    free(url);
    fclose(fp);

    if (urlList->head == NULL) {
        ListFree(urlList);
        return NULL;
    }

    return urlList;
}

Graph createGraph(List l) {
    Graph g = GraphNew(l->size);

    // Allocate memory the url filename string
    char *filename = malloc((MAX_STRLEN + strlen(".txt") + 1) * sizeof(char));
    if (filename == NULL) {
        fprintf(stderr, "error: out of memory\n");
        exit(EXIT_FAILURE);
    }
    
    char *url = malloc(MAX_STRLEN * sizeof(char) + 1);
    if (url == NULL) {
        fprintf(stderr, "error: out of memory\n");
        exit(EXIT_FAILURE);
    }

    for (Node curr = l->head; curr != NULL; curr = curr->next) {
        // Create a url node for each filename
        strcpy(filename, curr->url);
        strcat(filename, ".txt");

        // Open the url's file
        FILE *fp = fopen(filename, "r");
        if (fp == NULL) {
            fprintf(stderr, "fopen");
            exit(EXIT_FAILURE);
        }
        // Insert all the edges
        insertEdges(g, l, fp, url, curr);

    }
    free(filename);
    free(url);
    return g;
}

// Reads a url.txt file and inserts all edges the url is adjacent to
static void insertEdges(Graph g, List l, FILE *fp, char *url, Node curr) {
    // Move fp to be after the first line
    char buffer[MAX_STRLEN];
    fgets(buffer, MAX_STRLEN, fp);
    int file_offset = strlen(buffer);
    fseek(fp, file_offset, SEEK_SET);

    // Read strings until we read #end
    while (fscanf(fp, "%s", url) != EOF) {
        if (strcmp(url, "#end") == 0) break;
        
        int outlinkIndex = getUrlIndex(l, url);
        if (outlinkIndex >= 0 && outlinkIndex != curr->index) {
            Edge e = {curr->index, outlinkIndex, 1};
            GraphInsertEdge(g, e);
        }
    }
}

// Returns the index for url string
static int getUrlIndex(List l, char *url) {
    for (Node n = l->head; n != NULL; n = n->next) {
        if (strcmp(url, n->url) == 0) {
            return n->index;
        }
    }
    return -1;
}

// Calculates page ranks for each url using the given formula
List calculatePageRank(List l, Graph g, double d, double diffPR, 
                       int maxIterations) {
    double N = g->nV;
    
    // calculate iteration 0 rank, incoming and outgoing degree
    initialiseRankAndDegree(l, g, N);

    // Create graphs containing the Win and Wout values for each edge
    Graph gWin = setGraphWin(l, g);
    Graph gWout = setGraphWout(l, g);

    double diff = diffPR;
    for (int i = 1; i < maxIterations && diff >= diffPR; i++) {
        // store the previous rank
        for (Node pi = l->head; pi != NULL; pi = pi->next) {
            pi->prevRank = pi->rank;
        }
        // Update the rank
        for (Node pi = l->head; pi != NULL; pi = pi->next) {

            double weights =  getPageWeight(l, g, gWin, gWout, pi);
            pi->rank = (1 - d) / N + d * weights;
        }
        diff = calculateDiff(l);
    }
    GraphFree(gWin);
    GraphFree(gWout);
    return l;
}

// Initialises rank to 1/N and sets out/in degree for each url
static void initialiseRankAndDegree(List l, Graph g, double N) {
    for (Node pi = l->head; pi != NULL; pi = pi->next) {
        pi->rank = 1 / N;
        pi->outDegree = outgoingDegree(g, pi->index);
        pi->inDegree = incomingDegree(g, pi->index);
    }
}

// Takes a url index and returns the corresponding url node
static Node getNode(List l, int index) {
    assert(index >= 0 && index < l->size); 

    for (Node n = l->head; n != NULL; n = n->next) {
        if (n->index == index) return n;
    }
    return NULL;
}

// Sets the Win for each edge in the graph
static Graph setGraphWin(List l, Graph g) {
    Graph gWin = GraphNew(g->nV);
    for (int pj = 0; pj < g->nV; pj++) {
        for (int pi = 0; pi < g->nV; pi++) {
            // don't include self loops
            if (isAdjacent(g, pj, pi) && pj != pi) {
                Node pjNode = getNode(l, pj);
                Node piNode = getNode(l, pi);
                double Win = calculateWin(l, g, pjNode, piNode);

                Edge e = {pj, pi, Win};
                GraphInsertEdge(gWin, e);
            }
        }
    }
    return gWin;
}

// Sets the Wout for each edge in the graph
static Graph setGraphWout(List l, Graph g) {
    Graph gWout = GraphNew(g->nV);
    for (int pj = 0; pj < g->nV; pj++) {
        for (int pi = 0; pi < g->nV; pi++) {
            if (isAdjacent(g, pj, pi) && pj != pi) {
                Node pjNode = getNode(l, pj);
                Node piNode = getNode(l, pi);
                double Wout = calculateWout(l, g, pjNode, piNode);

                Edge e = {pj, pi, Wout};
                GraphInsertEdge(gWout, e);
            }
        }
    }
    return gWout;
}

// Calculates the Win for each edge
static double calculateWin(List l, Graph g, Node pj, Node pi) {
    // Store all the incoming links 
    double pjTotalIncomingLinks = 0;
    int pjRow = pj->index;
    for (int pjCol = 0; pjCol < g->nV; pjCol++) {
        if (isAdjacent(g, pjRow, pjCol)) {
            Node refPage = getNode(l, pjCol);
            pjTotalIncomingLinks += refPage->inDegree;
        }
    }
    return pi->inDegree / pjTotalIncomingLinks;
}

// Returns true if edges[v][w] == 1, else false.
static bool isAdjacent(Graph g, int v, int w) {
    if (g->edges[v][w] != 0) return true;
    return false; 
}

// Returns the number of incoming links for a url
static double incomingDegree(Graph g, int urlA) {
    int degree = 0;
    for (int urlB = 0; urlB < g->nV; urlB++) {
        // for each urlB check if it is an incoming link in urlA (i.e. urlB -> urlA)
        if (isAdjacent(g, urlB, urlA)) {
            degree++;
        }
    }
    return degree;
}

// Returns the number of outgoing links for a url
static double outgoingDegree(Graph g, int urlA) {
    double degree = 0;

    for (int urlB = 0; urlB < g->nV; urlB++) {
        if (isAdjacent(g, urlA, urlB)) {
            degree++;
        }
    }
    return degree;
}

// Calculates the Wout for an edge
static double calculateWout(List l, Graph g, Node pj, Node pi) {

    // Set to 0.5 if pi == 0 per spec
    double piOutgoingLinks = pi->outDegree;
    if (piOutgoingLinks == 0) {
        piOutgoingLinks = 0.5;
    }

    double pjTotalOutgoingLinks = 0;
    int pjRow = pj->index;
    for (int pjCol = 0; pjCol < g->nV; pjCol++) {

        if (isAdjacent(g, pjRow, pjCol)) {
            Node refPage = getNode(l, pjCol);
            if (refPage->outDegree == 0) {
                pjTotalOutgoingLinks += 0.5;
            } else {
                pjTotalOutgoingLinks += refPage->outDegree;
            }
        }
    }
    // Set to 0.5 per spec if == 0
    if (pjTotalOutgoingLinks == 0) {
        pjTotalOutgoingLinks = 0.5;
    } 
    return piOutgoingLinks / pjTotalOutgoingLinks;
}

// Calculates the page weight 
static double getPageWeight(List l, Graph g, Graph gWin, Graph gWout, Node pi) {
    int piIndex = pi->index;
    double totalWeight = 0;
    for (int pjIndex = 0; pjIndex < g->nV; pjIndex++) {
        if (isAdjacent(g, pjIndex, piIndex)) {
            Node pj = getNode(l, pjIndex);
            double Wout = gWout->edges[pjIndex][piIndex];
            double Win = gWin->edges[pjIndex][piIndex];
            double weight = pj->prevRank * Wout * Win;
            totalWeight += weight;
        }
    }
    return totalWeight;
}

// Calculates the diff for each iteration cycle
static double calculateDiff(List l) {
    double diff = 0;
    for (Node pi = l->head; pi != NULL; pi = pi->next) {
        diff += fabs(pi->rank - pi->prevRank);
    }
    return diff;
}