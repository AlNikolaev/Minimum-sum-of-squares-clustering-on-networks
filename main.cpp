#include <fstream>
#include <stdio.h>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <random>
#include <time.h>
#include <iostream>

using namespace std;

struct Edge {
	bool tabu; // true - if edge is tabu
	int firstEndpointIndex;
	int secondEndpointIndex;
	double length;
};

struct Interior {
	int edgeIndex;	
	double length; // length  belongs to (0;length of current edge); length - distance between firstEndpoint of edge and interior point 
};

struct Vertex {
	bool tabu; // true - if vertex is tabu
    int clusterIndex;
    double distanceBetweenVertexAndPrototype;
	double bottleneckPoint;
    double secondMinimumDistance;
	int clusterIndexOfSecondMinimumDistance;
};

struct IndexesArray {
	int * items;
	int size;
};

struct Prototype {
	double clusterObjective;
    // (vertexIndex = -1 and interior.edgeIndex != -1) or (vertexIndex != -1 and interior.edgeIndex = -1) 
	int vertexIndex;  
	Interior interior;
	IndexesArray assignedVertices;
};

struct ShortestPathToInteriorPoint {
	IndexesArray pathThroughFirstEndpoint;
	IndexesArray pathThroughSecondEndpoint;
	IndexesArray pathThroughDependOn_x;
	double sum_pathThroughFirstEndpoint;
	double sum_pathThroughSecondEndpoint;
};

struct BestOneClusterSolution {
	int edgeIndex;
	double e;
	double objective;
};

struct Clustering {
	Prototype* prototypeArray;
	double objective;
};

struct Interval {
	double endpoint;
	bool isLeftEndpoint;
	int vertex;
};


struct IntervalArray {
	Interval * items;
	int size;
};

const double MAX_DOUBLE = std::numeric_limits<double>::max();

int verticesNumber;
int edgesNumber;
int clustersNumber;
int counter;
int timeLimit;

Edge* edgeArray = 0;
Vertex* vertexArray = 0;
Clustering bestSolution;
Clustering currentSolution;
Clustering tmpSolution;
double ** shortestPathMatrix = 0;
int** isEdgeExist = 0;
ShortestPathToInteriorPoint shPath;  
double* segmentEnds = 0;
IndexesArray potentialAssignedVert;
IntervalArray intervalArray;
time_t startTime, endTime;
    
static int compare (const void* a, const void* b) {
	if (*(double *) b < *(double *) a) { return 1;}
	if (*(double *) b > *(double *) a) { return -1;}
	return 0;
}

static int compareInterval (const void* a, const void* b) {
	double aEndpoint = ((Interval *) a)->endpoint;
	double bEndpoint = ((Interval *) b)->endpoint;
	if (aEndpoint < bEndpoint) { return -1;}
	if (aEndpoint > bEndpoint) { return 1;}
	aEndpoint = ((Interval *) a)->vertex;
	bEndpoint = ((Interval *) b)->vertex;
	if (aEndpoint < bEndpoint) { return -1;}
	if (aEndpoint > bEndpoint) { return 1;}
	return 0;
}

void exitOnError(const char * error, const char * what)
{
	printf ("%s %s \n", error, what); 
	exit(1); 
}

/*find the lengths of the shortest paths between all pairs of vertices*/
void Floyd–Warshall_algorithm() 
{	
	for (int k = 1; k <= verticesNumber; k++) {
		for (int i = 1; i <= verticesNumber; i++) {
			for (int j = 1; j <= verticesNumber; j++) {
				if (shortestPathMatrix[i][k] != MAX_DOUBLE && shortestPathMatrix[k][j] != MAX_DOUBLE && shortestPathMatrix[i][j] > shortestPathMatrix[i][k] + shortestPathMatrix[k][j]) {
					shortestPathMatrix[i][j] = shortestPathMatrix[i][k] + shortestPathMatrix[k][j];
				}
			}
		}
	}
}

void readInstanceAndInitializeVariables(char name[]) 
{
	ifstream f (name);
	char buffer[256];
	if (! f.is_open()) 
	{ 
		exitOnError("Error opening file: ", name); 
	}
	
	do 
	{
		f.getline(buffer, 250);
	}
	while (buffer[0] == 'c' || buffer[0] == '\0'); 
	sscanf(buffer,"%d %d %d",&verticesNumber, &edgesNumber, &clustersNumber);

	if (verticesNumber < 1 || edgesNumber < 1 || clustersNumber < 1) {
		exitOnError("Wrong number of vertices or/and number of edges or/and number of clusters", ""); 
	}
	if (verticesNumber < clustersNumber) {
		exitOnError("Number of clusters should be less or equal to number of vertices", ""); 
	}
	edgeArray = new Edge[edgesNumber + 1];
	vertexArray = new Vertex[verticesNumber + 1];
	bestSolution.prototypeArray = new Prototype[clustersNumber + 1];
	bestSolution.objective = MAX_DOUBLE;
	currentSolution.prototypeArray = new Prototype[clustersNumber + 1];
	currentSolution.objective = MAX_DOUBLE;
	tmpSolution.prototypeArray = new Prototype[clustersNumber + 1];
	tmpSolution.objective = MAX_DOUBLE;
	shortestPathMatrix = new double* [verticesNumber + 1];
	isEdgeExist = new int * [verticesNumber + 1];
	segmentEnds = new double [2*verticesNumber];
	shPath.pathThroughDependOn_x.items = new int [verticesNumber];
	shPath.pathThroughFirstEndpoint.items = new int [verticesNumber];
	shPath.pathThroughSecondEndpoint.items = new int [verticesNumber];
    potentialAssignedVert.items = new int [verticesNumber];
	intervalArray.items = new Interval [2*verticesNumber];
	shPath.pathThroughDependOn_x.size = 0;
	shPath.pathThroughFirstEndpoint.size = 0;
	shPath.pathThroughSecondEndpoint.size = 0;
	shPath.sum_pathThroughFirstEndpoint = 0;
	shPath.sum_pathThroughSecondEndpoint = 0;
    potentialAssignedVert.size =  0;
	intervalArray.size = 0;
    for (int i = 1; i <= verticesNumber; i++) {
		shortestPathMatrix[i] = new double [verticesNumber + 1];
		vertexArray[i].tabu = false;
		vertexArray[i].distanceBetweenVertexAndPrototype = MAX_DOUBLE;
        isEdgeExist[i] = new int [verticesNumber + 1];
		for (int j = 1; j <= verticesNumber; j++) {
			isEdgeExist[i][j] = -1;
		}
		if (i <= clustersNumber) {
			bestSolution.prototypeArray[i].interior.edgeIndex = -1;
			bestSolution.prototypeArray[i].interior.length = 0;
			bestSolution.prototypeArray[i].vertexIndex = -1;
            bestSolution.prototypeArray[i].assignedVertices.items = new int [verticesNumber];
            bestSolution.prototypeArray[i].assignedVertices.size = 0;
			
			currentSolution.prototypeArray[i].interior.edgeIndex = -1;
			currentSolution.prototypeArray[i].interior.length = 0;
			currentSolution.prototypeArray[i].vertexIndex = -1;
            currentSolution.prototypeArray[i].assignedVertices.items = new int [verticesNumber];
            currentSolution.prototypeArray[i].assignedVertices.size = 0;
			
			tmpSolution.prototypeArray[i].interior.edgeIndex = -1;
			tmpSolution.prototypeArray[i].interior.length = 0;
			tmpSolution.prototypeArray[i].vertexIndex = -1;
            tmpSolution.prototypeArray[i].assignedVertices.items = new int [verticesNumber];
            tmpSolution.prototypeArray[i].assignedVertices.size = 0;
			
			for (int j = 0; j < verticesNumber; j++) {
				 bestSolution.prototypeArray[i].assignedVertices.items[j] = -1;
				 currentSolution.prototypeArray[i].assignedVertices.items[j] = -1;
			}
		}
		for (int j = 1; j <= verticesNumber; j++) {
			if (i != j) {
				shortestPathMatrix[i][j] = MAX_DOUBLE;
			}
			else {
				shortestPathMatrix[i][j] = 0;
			}
		}
	}

	int edgeIndex = 1;
	while (! f.eof())
	{
		f.getline(buffer, 250);
		if (buffer[0] == 'c' || buffer[0] == '\0') continue;
		int result = sscanf(buffer,"%d %d %lf", &edgeArray[edgeIndex].firstEndpointIndex, &edgeArray[edgeIndex].secondEndpointIndex, &edgeArray[edgeIndex].length);
		if (result != 3) {
			exitOnError("Wrong line: ", buffer); 
		}
		if (isEdgeExist[edgeArray[edgeIndex].firstEndpointIndex][edgeArray[edgeIndex].secondEndpointIndex] == -1) {
			isEdgeExist[edgeArray[edgeIndex].firstEndpointIndex][edgeArray[edgeIndex].secondEndpointIndex] = edgeIndex;
			isEdgeExist[edgeArray[edgeIndex].secondEndpointIndex][edgeArray[edgeIndex].firstEndpointIndex] = edgeIndex;
			if (edgeArray[edgeIndex].firstEndpointIndex > edgeArray[edgeIndex].secondEndpointIndex) {
				int tmp = edgeArray[edgeIndex].firstEndpointIndex;
				edgeArray[edgeIndex].firstEndpointIndex = edgeArray[edgeIndex].secondEndpointIndex;
				edgeArray[edgeIndex].secondEndpointIndex = tmp;
			}
			edgeArray[edgeIndex].tabu = false;
			shortestPathMatrix[edgeArray[edgeIndex].firstEndpointIndex][edgeArray[edgeIndex].secondEndpointIndex] = edgeArray[edgeIndex].length;
			shortestPathMatrix[edgeArray[edgeIndex].secondEndpointIndex][edgeArray[edgeIndex].firstEndpointIndex] = edgeArray[edgeIndex].length;
			edgeIndex++;
		}
		else {
			int curEdgeIndex = isEdgeExist[edgeArray[edgeIndex].firstEndpointIndex][edgeArray[edgeIndex].secondEndpointIndex];
			if (edgeArray[edgeIndex].length < edgeArray[curEdgeIndex].length) {
				shortestPathMatrix[edgeArray[curEdgeIndex].firstEndpointIndex][edgeArray[curEdgeIndex].secondEndpointIndex] = edgeArray[edgeIndex].length;
				shortestPathMatrix[edgeArray[curEdgeIndex].secondEndpointIndex][edgeArray[curEdgeIndex].firstEndpointIndex] = edgeArray[edgeIndex].length;
			}
			edgesNumber--;
		}
	}

	if (edgeIndex != (edgesNumber + 1)) {
		exitOnError("Wrong number of edges", ""); 
	}
	double density = (double) edgesNumber / (verticesNumber * (verticesNumber - 1) / 2);
	printf("|E| = %d |V| = %d density = %6.5f  \n", edgesNumber, verticesNumber, density);
    f.close();
}

std::random_device rd;
std::mt19937_64 generator(rd());
std::uniform_int_distribution<int>* vertexDiscreteDistribution = NULL;
std::uniform_int_distribution<int>* edgeDiscreteDistribution = NULL;
std::uniform_real_distribution<double> continuousDistributionWithout0and1(0.00000000000000000000000000000000001,1);
std::uniform_real_distribution<double> continuousDistributionWith0and1(0,1.00000000000000000000000000000000001);
std::uniform_int_distribution<int> coin(0, 1);

void generate_prototype_vertex_with_max_dist(int clusterIndex, Prototype* prototypeArray);

void generate_prototype_feasible(int clusterIndex, Prototype* prototypeArray) { 
	int currentEdge;
    counter = 0;
	do {
		currentEdge = (*edgeDiscreteDistribution)(generator);
        counter++;
    }
    while ( (edgeArray[currentEdge].tabu || (vertexArray[edgeArray[currentEdge].firstEndpointIndex].tabu && vertexArray[edgeArray[currentEdge].secondEndpointIndex].tabu)) && counter < 1000 );

    if (counter > 1000) {
        cout << "generate_prototype_vertex_with_max_dist from generate_prototype_feasible!\n";
        generate_prototype_vertex_with_max_dist(clusterIndex, prototypeArray);
    }
    else {
        double x = continuousDistributionWith0and1(generator);
        if ( (!vertexArray[edgeArray[currentEdge].firstEndpointIndex].tabu && vertexArray[edgeArray[currentEdge].secondEndpointIndex].tabu) || x == 0) {
            vertexArray[edgeArray[currentEdge].firstEndpointIndex].tabu = true;
            prototypeArray[clusterIndex].vertexIndex = edgeArray[currentEdge].firstEndpointIndex;
            prototypeArray[clusterIndex].interior.edgeIndex = -1;
        }
        else if ( (!vertexArray[edgeArray[currentEdge].secondEndpointIndex].tabu && vertexArray[edgeArray[currentEdge].firstEndpointIndex].tabu) || x >= 1) {
            vertexArray[edgeArray[currentEdge].secondEndpointIndex].tabu = true;
            prototypeArray[clusterIndex].vertexIndex = edgeArray[currentEdge].secondEndpointIndex;
            prototypeArray[clusterIndex].interior.edgeIndex = -1;
        } else {
            prototypeArray[clusterIndex].interior.edgeIndex = currentEdge;
            prototypeArray[clusterIndex].interior.length = x * edgeArray[currentEdge].length;
            prototypeArray[clusterIndex].vertexIndex = -1;
            edgeArray[currentEdge].tabu = true;
            vertexArray[edgeArray[currentEdge].firstEndpointIndex].tabu = true;
            vertexArray[edgeArray[currentEdge].secondEndpointIndex].tabu = true;
        }    
    }
}

void generate_prototype_inside_edge(int clusterIndex, Prototype* prototypeArray) {
	int currentEdge;
    do {
		currentEdge = (*edgeDiscreteDistribution)(generator);
    }
    while ( edgeArray[currentEdge].tabu || vertexArray[edgeArray[currentEdge].firstEndpointIndex].tabu || vertexArray[edgeArray[currentEdge].secondEndpointIndex].tabu );

    double x = continuousDistributionWithout0and1(generator);
	prototypeArray[clusterIndex].interior.edgeIndex = currentEdge;
	prototypeArray[clusterIndex].interior.length = x * edgeArray[currentEdge].length;
	prototypeArray[clusterIndex].vertexIndex = -1;
	edgeArray[currentEdge].tabu = true;
    vertexArray[edgeArray[currentEdge].firstEndpointIndex].tabu = true;
    vertexArray[edgeArray[currentEdge].secondEndpointIndex].tabu = true;
}

void generate_prototype_on_vertex(int clusterIndex, Prototype* prototypeArray) {
	int currentVertex;
    do {
		currentVertex = (*vertexDiscreteDistribution)(generator);
    }
    while ( vertexArray[currentVertex].tabu );

	vertexArray[currentVertex].tabu = true;
	prototypeArray[clusterIndex].vertexIndex = currentVertex;
	prototypeArray[clusterIndex].interior.edgeIndex = -1;
}

void generate_prototype_on_vertex_or_edge_50_50(int clusterIndex, Prototype* prototypeArray) {
	if (coin(generator)) {
		generate_prototype_on_vertex(clusterIndex, prototypeArray);
	}
    else {
		generate_prototype_inside_edge(clusterIndex, prototypeArray);
	}    
}

void generate_prototype_on_sorted_edges(Prototype* prototypeArray) { // it is possible to make it more effective
	int* indexes = new int [edgesNumber + 1];  
	int i, j;
	for (i = 1; i <= edgesNumber; i++) {
		indexes[i] = i;
	}
	for (i = 1; i <= edgesNumber; i++) { // sort indexes by increasing edge length
		for (j = i; j <= edgesNumber; j++) {
			if (edgeArray[indexes[i]].length > edgeArray[indexes[j]].length) {
				int tmp;
				tmp = indexes[j];
				indexes[j] = indexes[i];
				indexes[i] = tmp;
			}			
		}	
	}
	int startIndex = 1;
	int currentEdge = 0;
	for (i = 1; i <= clustersNumber; i++) {
		int counter = 0;
		for (j = startIndex; (edgeArray[indexes[j]].length == edgeArray[indexes[startIndex]].length) && j <= edgesNumber; j++) 
		{
			if (!edgeArray[indexes[j]].tabu && !vertexArray[edgeArray[indexes[j]].firstEndpointIndex].tabu && !vertexArray[edgeArray[indexes[j]].secondEndpointIndex].tabu) {
				counter++;
			}
		}
		if (counter == 0) {
			startIndex = j;
			i--;
			continue;
		}
		if (startIndex != (j - 1)) {
			std::uniform_int_distribution<int> discreteDistribution(startIndex, j - 1);
			do {
				currentEdge = indexes [discreteDistribution(generator)];
			}
			while ( edgeArray[currentEdge].tabu || vertexArray[edgeArray[currentEdge].firstEndpointIndex].tabu || vertexArray[edgeArray[currentEdge].secondEndpointIndex].tabu );
		}
		else {
			currentEdge = indexes [startIndex];
		}
		
		double x = continuousDistributionWithout0and1(generator);
		prototypeArray[i].interior.edgeIndex = currentEdge;
		prototypeArray[i].interior.length = x * edgeArray[currentEdge].length;
		prototypeArray[i].vertexIndex = -1;
		edgeArray[currentEdge].tabu = true;
		vertexArray[edgeArray[currentEdge].firstEndpointIndex].tabu = true;
		vertexArray[edgeArray[currentEdge].secondEndpointIndex].tabu = true;
		if (counter == 1) {
			startIndex = j;
		}
	}
}

void generate_prototype(char generatorName, int clusterIndex, Prototype* prototypeArray) {
	switch(generatorName) {
		case 'v': {
			generate_prototype_on_vertex(clusterIndex, prototypeArray);
			break;
		}
		case 'e': {
			generate_prototype_inside_edge(clusterIndex, prototypeArray);
			break;
		}
		case 'f': {
			generate_prototype_feasible(clusterIndex, prototypeArray);
			break;		  
		}
		case 'r': {
			generate_prototype_on_vertex_or_edge_50_50(clusterIndex, prototypeArray);
			break;		 
		}
        default: {
            exitOnError("Wrong generator parameter: ", &generatorName);         
        }
	}
}

void generate_prototype_vertex_with_max_dist(int clusterIndex, Prototype* prototypeArray) {
    int indexOfVertexWithMaxDist = 1;
    for (int vertexIndex = 2; vertexIndex <= verticesNumber; vertexIndex++) {
        if (vertexArray[vertexIndex].distanceBetweenVertexAndPrototype > vertexArray[indexOfVertexWithMaxDist].distanceBetweenVertexAndPrototype) {
            indexOfVertexWithMaxDist = vertexIndex;
        }
    }
    prototypeArray[clusterIndex].interior.edgeIndex = -1;
    prototypeArray[clusterIndex].vertexIndex = indexOfVertexWithMaxDist;
    vertexArray[indexOfVertexWithMaxDist].distanceBetweenVertexAndPrototype = 0;
}

void generate_prototype_for_remove_degeneracy(char generatorName, int clusterIndex, Prototype* prototypeArray) {
	if (generatorName == 'm') {
		generate_prototype_vertex_with_max_dist(clusterIndex, prototypeArray);
	}
	else {
        generate_prototype(generatorName, clusterIndex, prototypeArray);
	}
}

void generate_initialSolution(char generatorName, Prototype* prototypeArray) {
	if (generatorName == 's') {
		generate_prototype_on_sorted_edges(prototypeArray);
	}
	else {
		for (int i = 1; i <= clustersNumber; i++) {
			generate_prototype(generatorName, i, prototypeArray);
		}
	}
}

/* return true if prototypes are changed*/
bool allocate(Clustering& currentSolution) {
	Prototype* prototypeArray = currentSolution.prototypeArray;
	currentSolution.objective = 0;
	bool functionResult = false;
	double minDistance;
    double currentDistance;
	for (int cluster = 1; cluster <= clustersNumber; cluster++) {
		prototypeArray[cluster].clusterObjective = 0;
        prototypeArray[cluster].assignedVertices.size = 0;
	}
    for (int currentVertex = 1; currentVertex <= verticesNumber; currentVertex++) {
        minDistance = MAX_DOUBLE;
        int nearestPrototypeIndex = -1;
        // find the nearest prototype
        for (int cluster = 1; cluster <= clustersNumber; cluster++) { 
            if (prototypeArray[cluster].vertexIndex != -1) {
                currentDistance = shortestPathMatrix[currentVertex][prototypeArray[cluster].vertexIndex];
            }
            else {
                int edgeIndex = prototypeArray[cluster].interior.edgeIndex;
				double lengthFromFirstendpoint = prototypeArray[cluster].interior.length;
                
                currentDistance = shortestPathMatrix[currentVertex][edgeArray[edgeIndex].firstEndpointIndex] + lengthFromFirstendpoint;
                double tmpDistance = shortestPathMatrix[currentVertex][edgeArray[edgeIndex].secondEndpointIndex] + edgeArray[edgeIndex].length - lengthFromFirstendpoint;
                if (tmpDistance < currentDistance) {
                    currentDistance = tmpDistance;
                }
            }
            if (minDistance > currentDistance) {
                minDistance = currentDistance;
                nearestPrototypeIndex = cluster;
            }
        }
		if (prototypeArray[nearestPrototypeIndex].assignedVertices.items[prototypeArray[nearestPrototypeIndex].assignedVertices.size] != currentVertex) 
		{
			functionResult = true;
		}
        prototypeArray[nearestPrototypeIndex].assignedVertices.items[prototypeArray[nearestPrototypeIndex].assignedVertices.size++] = currentVertex;
        vertexArray[currentVertex].clusterIndex = nearestPrototypeIndex;
        vertexArray[currentVertex].distanceBetweenVertexAndPrototype = minDistance;
        prototypeArray[nearestPrototypeIndex].clusterObjective += pow(minDistance, 2);
	}
    for (int cluster = 1; cluster <= clustersNumber; cluster++) { 
        prototypeArray[cluster].assignedVertices.items[prototypeArray[cluster].assignedVertices.size] = -1;
		currentSolution.objective += prototypeArray[cluster].clusterObjective;
	}
	return functionResult;
}

inline double computeDistance(Prototype& currentPrototype, int currentVertex) {
	if (currentPrototype.vertexIndex != -1) {
		return shortestPathMatrix[currentVertex][currentPrototype.vertexIndex];
	}
	else {
		int edgeIndex = currentPrototype.interior.edgeIndex;
		double lengthFromFirstendpoint = currentPrototype.interior.length;
		return min(shortestPathMatrix[currentVertex][edgeArray[edgeIndex].firstEndpointIndex] + lengthFromFirstendpoint, 
			shortestPathMatrix[currentVertex][edgeArray[edgeIndex].secondEndpointIndex] + edgeArray[edgeIndex].length - lengthFromFirstendpoint);
	}
}

inline void updateTabuInformation(Prototype& currentPrototype) {
	if (currentPrototype.vertexIndex != -1) {
		vertexArray[currentPrototype.vertexIndex].tabu = true;
	}
    else {
		int edgeIndex = currentPrototype.interior.edgeIndex;
        edgeArray[edgeIndex].tabu = true;
        vertexArray[edgeArray[edgeIndex].firstEndpointIndex].tabu = true;
        vertexArray[edgeArray[edgeIndex].secondEndpointIndex].tabu = true;
    }
}

inline void clearTabuInformation(Prototype& currentPrototype) {
	if (currentPrototype.vertexIndex != -1) {
		vertexArray[currentPrototype.vertexIndex].tabu = false;
	}
    else {
		int edgeIndex = currentPrototype.interior.edgeIndex;
        edgeArray[edgeIndex].tabu = false;
        vertexArray[edgeArray[edgeIndex].firstEndpointIndex].tabu = false;
        vertexArray[edgeArray[edgeIndex].secondEndpointIndex].tabu = false;
    }
}

inline void clearTabuInformationForNotRemaining(bool* remainingPrototypes, Prototype* prototypeArray) {
    for (int cluster = 1; cluster <= clustersNumber; cluster++) {
        if (!remainingPrototypes[cluster]) {
			clearTabuInformation(prototypeArray[cluster]);
		}
    }
}

inline void clearAllTabuInformation() {
    for (int i = 1; i <= verticesNumber; i++) {
        vertexArray[i].tabu = false;
    }
    for (int i = 1; i <= edgesNumber; i++) {
        edgeArray[i].tabu = false;
    }
}

void shaking(bool* remainingPrototypes, int numberOfDifferencesBetweenSolutions, Prototype* prototypeArray, char generatorName) {
    std::random_device rd;
    std::mt19937 generator(rd());
    std::uniform_int_distribution<int> discreteDistribution(1, clustersNumber);
    int currentNumberOfDiff;
    int currentCluster;
    if (numberOfDifferencesBetweenSolutions > (clustersNumber/2)) {
        for (int cluster = 1; cluster <= clustersNumber; cluster++) {
            remainingPrototypes[cluster] = false;
        }
        currentNumberOfDiff = clustersNumber; 
        while (currentNumberOfDiff > numberOfDifferencesBetweenSolutions) {
            currentCluster = discreteDistribution(generator);
            if (!remainingPrototypes[currentCluster]) {
                remainingPrototypes[currentCluster] = true;
                currentNumberOfDiff--;
            }
        }
    } else {
        for (int cluster = 1; cluster <= clustersNumber; cluster++) {
            remainingPrototypes[cluster] = true;
        }
        currentNumberOfDiff = 0; 
        while (currentNumberOfDiff < numberOfDifferencesBetweenSolutions) {
            currentCluster = discreteDistribution(generator);
            if (remainingPrototypes[currentCluster]) {
                remainingPrototypes[currentCluster] = false;
                currentNumberOfDiff++;
            }
        }

    }
    clearTabuInformationForNotRemaining(remainingPrototypes, prototypeArray);
    for (int cluster = 1; cluster <= clustersNumber; cluster++) {
        if (!remainingPrototypes[cluster]) {
			//generate prototype
			generate_prototype(generatorName, cluster, prototypeArray);
        }
    }
}

int remove_degeneracy(bool* remainingPrototypes, Prototype* prototypeArray, char generatorName) {
    int counterDegeneracy = 0;
    for (int cluster = 1; cluster <= clustersNumber; cluster++) {
        if (prototypeArray[cluster].assignedVertices.size == 0) {
            clearTabuInformation(prototypeArray[cluster]);
			++counterDegeneracy;
        }
    }
    if (counterDegeneracy == 0) { return 0;}
	for (int cluster = 1; cluster <= clustersNumber; cluster++) {
        if (prototypeArray[cluster].assignedVertices.size == 0) {
            remainingPrototypes[cluster] = false;
			//generate prototype
			generate_prototype_for_remove_degeneracy(generatorName, cluster, prototypeArray);
        }
    }
    return counterDegeneracy;
}

void printEdge(int index) {
    cout << "(" << edgeArray[index].firstEndpointIndex << "; " << edgeArray[index].secondEndpointIndex << ")";
}

void printPrototypes(Prototype* prototypeArray) {
    for (int i = 1; i <= clustersNumber; i++) {
        cout << "Prototype " << i << " ";
		if (prototypeArray[i].vertexIndex != -1) {
            cout << "vertex " << prototypeArray[i].vertexIndex;     
        }
        else {
            cout << "edge ";
			printEdge(prototypeArray[i].interior.edgeIndex);
			cout << " x=" << (prototypeArray[i].interior.length / edgeArray[prototypeArray[i].interior.edgeIndex].length);
        }
        cout << '\n';
    }
}

double computeObjectiveFunctionForClusterAccordingToShPath(Prototype& currentPrototype, int edgeIndex, double e) {
	double sum = 0;
    if (e > 0 && e < edgeArray[edgeIndex].length) {
        for (int i = 0; i < shPath.pathThroughFirstEndpoint.size; i++) {
		    sum += pow((shortestPathMatrix[shPath.pathThroughFirstEndpoint.items[i]][edgeArray[edgeIndex].firstEndpointIndex] + e), 2);
	    }
	    for (int i = 0; i < shPath.pathThroughDependOn_x.size; i++) {
		    sum += pow((shortestPathMatrix[shPath.pathThroughDependOn_x.items[i]][edgeArray[edgeIndex].firstEndpointIndex] + e), 2);
	    }
	    e = edgeArray[edgeIndex].length - e;
	    for (int i = 0; i < shPath.pathThroughSecondEndpoint.size; i++) {
		    sum += pow((shortestPathMatrix[shPath.pathThroughSecondEndpoint.items[i]][edgeArray[edgeIndex].secondEndpointIndex] + e), 2);
	    }
    }
    else {
        int vertexPr;
        if (e <= 0) {
            vertexPr = edgeArray[edgeIndex].firstEndpointIndex;
        } else {
            vertexPr = edgeArray[edgeIndex].secondEndpointIndex;
        }
        for (int index = 0; index < currentPrototype.assignedVertices.size; index++) {
            sum += pow (shortestPathMatrix[vertexPr][currentPrototype.assignedVertices.items[index]], 2);
        }
    }
    return sum;
}

double computeObjectiveForClusterVertices(IndexesArray& assignedVertices, int edgeIndex, double e) {
	double sum = 0;
    double dist = 0;
    int vertexIndex;
    if (e > 0 && e < edgeArray[edgeIndex].length) {
        for (int index = 0; index < assignedVertices.size; index++) {
            vertexIndex = assignedVertices.items[index];
            dist = min(shortestPathMatrix[vertexIndex][edgeArray[edgeIndex].firstEndpointIndex] + e, shortestPathMatrix[vertexIndex][edgeArray[edgeIndex].secondEndpointIndex] + edgeArray[edgeIndex].length - e);
	        sum += pow (dist, 2);
        }
    }
    else {
        if (e <= 0) {
            vertexIndex = edgeArray[edgeIndex].firstEndpointIndex;
        } else {
            vertexIndex = edgeArray[edgeIndex].secondEndpointIndex;
        }
        for (int index = 0; index < assignedVertices.size; index++) {
            sum += pow (shortestPathMatrix[vertexIndex][assignedVertices.items[index]], 2);
        }
    }
    return sum;
}

double computeObjectiveForCluster(Prototype& currentPrototype) {
	double sum = 0;
    double dist = 0;
    int vertexIndex;
	if (currentPrototype.interior.edgeIndex != -1) {
		for (int index = 0; index < currentPrototype.assignedVertices.size; index++) {
            vertexIndex = currentPrototype.assignedVertices.items[index];
			dist = min(shortestPathMatrix[vertexIndex][edgeArray[currentPrototype.interior.edgeIndex].firstEndpointIndex] + currentPrototype.interior.length, shortestPathMatrix[vertexIndex][edgeArray[currentPrototype.interior.edgeIndex].secondEndpointIndex] + edgeArray[currentPrototype.interior.edgeIndex].length - currentPrototype.interior.length);
	        sum += pow (dist, 2);
        }
    }
    else {
        for (int index = 0; index < currentPrototype.assignedVertices.size; index++) {
			sum += pow (shortestPathMatrix[currentPrototype.vertexIndex][currentPrototype.assignedVertices.items[index]], 2);
        }
    }
    return sum;
}

void computeBottleneckPointsAndSegmentEnds(IndexesArray assignedVertices, int currentEdge, double intervalPoint1, double intervalPoint2, int& sizeOfSegmentEndsArray) {
	shPath.pathThroughDependOn_x.size = 0;
	shPath.pathThroughFirstEndpoint.size = 0;
	shPath.pathThroughSecondEndpoint.size = 0;
	shPath.sum_pathThroughFirstEndpoint = 0;
	shPath.sum_pathThroughSecondEndpoint = 0;
    for (int currentIndex = 0; currentIndex < assignedVertices.size; currentIndex++) {
		int currentVertex = assignedVertices.items[currentIndex];
		vertexArray[currentVertex].bottleneckPoint = 0.5*(edgeArray[currentEdge].length + shortestPathMatrix[currentVertex][edgeArray[currentEdge].secondEndpointIndex] - shortestPathMatrix[currentVertex][edgeArray[currentEdge].firstEndpointIndex]);
		if (vertexArray[currentVertex].bottleneckPoint <= intervalPoint1) {
			shPath.pathThroughSecondEndpoint.items[shPath.pathThroughSecondEndpoint.size++] = currentVertex;
			shPath.sum_pathThroughSecondEndpoint += shortestPathMatrix[currentVertex][edgeArray[currentEdge].secondEndpointIndex];
		}
		else if (vertexArray[currentVertex].bottleneckPoint >= intervalPoint2) {
			shPath.pathThroughFirstEndpoint.items[shPath.pathThroughFirstEndpoint.size++] = currentVertex;
			shPath.sum_pathThroughFirstEndpoint += shortestPathMatrix[currentVertex][edgeArray[currentEdge].firstEndpointIndex];
		}
		else {
			shPath.pathThroughDependOn_x.items[shPath.pathThroughDependOn_x.size++] = currentVertex;
			segmentEnds[sizeOfSegmentEndsArray++] = vertexArray[currentVertex].bottleneckPoint; 
			shPath.sum_pathThroughFirstEndpoint += shortestPathMatrix[currentVertex][edgeArray[currentEdge].firstEndpointIndex];
		}
	}
}

bool solve_1_prototype_for_edge(int currentEdge, Prototype& currentPrototype, BestOneClusterSolution& bestSolution, double intervalPoint1,  double intervalPoint2) {
	if (intervalPoint1 > intervalPoint2 || intervalPoint1 > edgeArray[currentEdge].length || intervalPoint2 < 0) {
		return true;
	} else if (intervalPoint1 == intervalPoint2) {
		double result = computeObjectiveForClusterVertices(currentPrototype.assignedVertices, currentEdge, intervalPoint1);
		if (result < bestSolution.objective) {
			bestSolution.edgeIndex = currentEdge;
			bestSolution.e = intervalPoint1;
			bestSolution.objective = result;
		}
		return true;
	}

    // compute bottleneck points
    segmentEnds[0] = intervalPoint1;
	segmentEnds[1] = intervalPoint2;
	int sizeOfSegmentEndsArray = 2;
    computeBottleneckPointsAndSegmentEnds(currentPrototype.assignedVertices, currentEdge, intervalPoint1, intervalPoint2, sizeOfSegmentEndsArray);

    // sort segmentEnds by incresing order
	qsort(segmentEnds,sizeOfSegmentEndsArray, sizeof(double), compare);
	double e = 0;
	// find all local minimums
	for (int segmentIndex = 0; segmentIndex < (sizeOfSegmentEndsArray - 1); segmentIndex++) {
		if (segmentEnds[segmentIndex] == segmentEnds[segmentIndex + 1]) { continue; }
		if (segmentIndex != 0) {
			for (int index = 0; index < shPath.pathThroughDependOn_x.size; index++) {
				int currentVertex = shPath.pathThroughDependOn_x.items[index];
				if ( vertexArray[currentVertex].bottleneckPoint <= segmentEnds[segmentIndex]) {
					shPath.pathThroughSecondEndpoint.items[shPath.pathThroughSecondEndpoint.size++] = currentVertex;
					if (shPath.pathThroughDependOn_x.size > 1) {
						shPath.pathThroughDependOn_x.items[index] = shPath.pathThroughDependOn_x.items[shPath.pathThroughDependOn_x.size - 1]; 
					}
					--shPath.pathThroughDependOn_x.size; 
					shPath.sum_pathThroughFirstEndpoint -= shortestPathMatrix[currentVertex][edgeArray[currentEdge].firstEndpointIndex];
					shPath.sum_pathThroughSecondEndpoint += shortestPathMatrix[currentVertex][edgeArray[currentEdge].secondEndpointIndex];
					index--;
				}
			}	
		}
		double leftSum = shPath.sum_pathThroughFirstEndpoint + segmentEnds[segmentIndex] * (shPath.pathThroughFirstEndpoint.size + shPath.pathThroughDependOn_x.size);
		double rightSum = shPath.sum_pathThroughSecondEndpoint + (edgeArray[currentEdge].length - segmentEnds[segmentIndex]) * (shPath.pathThroughSecondEndpoint.size);
		if ( leftSum < rightSum) {
			if (segmentEnds[segmentIndex] + (rightSum - leftSum) / (double) currentPrototype.assignedVertices.size > segmentEnds[segmentIndex + 1]) {
				if ( (segmentIndex + 2) != sizeOfSegmentEndsArray ) {
					continue;
				}
				e = segmentEnds[segmentIndex + 1];
			} else { e = segmentEnds[segmentIndex] + (rightSum - leftSum) / (double) currentPrototype.assignedVertices.size;}
		} 
		else { 
			if (segmentIndex != 0) {continue;} 
			e = 0;
		}
		double result = computeObjectiveFunctionForClusterAccordingToShPath(currentPrototype, currentEdge, e);
        if (result < bestSolution.objective) {
			bestSolution.edgeIndex = currentEdge;
			bestSolution.e = e;
			bestSolution.objective = result;
		}
	}
	return true;
}

// change prototype
void changePrototypeToBestOneClusterSolution (Prototype& currentPrototype, BestOneClusterSolution& bestSolution) {
	currentPrototype.clusterObjective = bestSolution.objective;
    if (bestSolution.e > 0 && bestSolution.e < edgeArray[bestSolution.edgeIndex].length) {
		currentPrototype.vertexIndex = -1;
		currentPrototype.interior.edgeIndex = bestSolution.edgeIndex;
		currentPrototype.interior.length = bestSolution.e;
	}
	else if (bestSolution.e == 0) {
		currentPrototype.vertexIndex = edgeArray[bestSolution.edgeIndex].firstEndpointIndex;
		currentPrototype.interior.edgeIndex = -1;
	}
	else {
		currentPrototype.vertexIndex = edgeArray[bestSolution.edgeIndex].secondEndpointIndex;
		currentPrototype.interior.edgeIndex = -1;
	}
}

void solve_1_prototype_for_every_cluster(Clustering& currentSolution) {
	BestOneClusterSolution bestSolution;
    Prototype* prototypeArray = currentSolution.prototypeArray;
    currentSolution.objective = 0;
    for (int cluster = 1; cluster <= clustersNumber; cluster++) {
		clearTabuInformation(prototypeArray[cluster]);
		if (prototypeArray[cluster].assignedVertices.size > 1) {
			bestSolution.objective = prototypeArray[cluster].clusterObjective;
			int clusterSize = prototypeArray[cluster].assignedVertices.size;
			for (int i = 0; i < clusterSize; i++) {
				for (int j = 0; j < clusterSize; j++) {
					int currentEdge = isEdgeExist[prototypeArray[cluster].assignedVertices.items[i]][prototypeArray[cluster].assignedVertices.items[j]];
					if (currentEdge != -1) {
						solve_1_prototype_for_edge(currentEdge, prototypeArray[cluster], bestSolution, 0.0, edgeArray[currentEdge].length);
					}
				}
			}
			if (bestSolution.objective < prototypeArray[cluster].clusterObjective) {
				changePrototypeToBestOneClusterSolution(prototypeArray[cluster], bestSolution);
            }
		}
		else if (prototypeArray[cluster].assignedVertices.size == 1) {
			prototypeArray[cluster].vertexIndex = prototypeArray[cluster].assignedVertices.items[0];
			prototypeArray[cluster].interior.edgeIndex = -1;
			prototypeArray[cluster].clusterObjective = 0;
		}
		else {
			cout << "Error! solve_1_prototype_for_every_cluster (852 line)\n";
		}
		updateTabuInformation(prototypeArray[cluster]);
		double tmpResult = computeObjectiveForCluster(prototypeArray[cluster]); 
		if (abs (tmpResult - prototypeArray[cluster].clusterObjective) > 0.0001) {
			cout << "Error! solve_1_prototype_for_every_cluster (857 line)\n";
		}
		currentSolution.objective += prototypeArray[cluster].clusterObjective;
	}
}

/*solve 1-opt for cluster with deleted vertex*/
void solve_1_prototype_for_cluster_with_deleted_vertex(Prototype& currentPrototype, int deletedVertexIndex) {
	double dist = 0;
	if ( currentPrototype.interior.edgeIndex != -1) {
		dist = min(shortestPathMatrix[deletedVertexIndex][edgeArray[currentPrototype.interior.edgeIndex].firstEndpointIndex] + currentPrototype.interior.length, shortestPathMatrix[deletedVertexIndex][edgeArray[currentPrototype.interior.edgeIndex].secondEndpointIndex] + edgeArray[currentPrototype.interior.edgeIndex].length - currentPrototype.interior.length);
	} else {
		dist = shortestPathMatrix[deletedVertexIndex][currentPrototype.vertexIndex];
	}
    BestOneClusterSolution bestSolution;
    if (currentPrototype.assignedVertices.size > 1) {
		bestSolution.objective = currentPrototype.clusterObjective - dist * dist;
		currentPrototype.clusterObjective = bestSolution.objective;
		int clusterSize = currentPrototype.assignedVertices.size;
		for (int i = 0; i < clusterSize; i++) {
			for (int j = 0; j < clusterSize; j++) {
				int currentEdge = isEdgeExist[currentPrototype.assignedVertices.items[i]][currentPrototype.assignedVertices.items[j]];
				if (currentEdge != -1) {
					double bottleneckPoint = 0.5*(edgeArray[currentEdge].length + shortestPathMatrix[deletedVertexIndex][edgeArray[currentEdge].secondEndpointIndex] - shortestPathMatrix[deletedVertexIndex][edgeArray[currentEdge].firstEndpointIndex]);
					double epsilon = shortestPathMatrix[deletedVertexIndex][edgeArray[currentEdge].firstEndpointIndex] + abs(bottleneckPoint) - dist;
					double intervalPoint1, intervalPoint2; 
					if (bottleneckPoint > 0 && bottleneckPoint < edgeArray[currentEdge].length) {
						intervalPoint1 = max(0.0, bottleneckPoint - epsilon);
						intervalPoint2 = min(edgeArray[currentEdge].length, bottleneckPoint + epsilon);
					}
					else if (bottleneckPoint <= 0) {
						intervalPoint1 = max(0.0, bottleneckPoint + epsilon);
						intervalPoint2 = min(edgeArray[currentEdge].length, bottleneckPoint + epsilon);
					}
					else {
						intervalPoint1 = max(0.0, bottleneckPoint - epsilon);
						intervalPoint2 = min(edgeArray[currentEdge].length, bottleneckPoint - epsilon);
					}
					solve_1_prototype_for_edge(currentEdge, currentPrototype, bestSolution, intervalPoint1, intervalPoint2);
				}	
			}
		}
		if (bestSolution.objective < currentPrototype.clusterObjective) {
        	changePrototypeToBestOneClusterSolution(currentPrototype, bestSolution);
        }
	} 
	else if (currentPrototype.assignedVertices.size == 1) {
		currentPrototype.vertexIndex = currentPrototype.assignedVertices.items[0];
		currentPrototype.interior.edgeIndex = -1;
		currentPrototype.clusterObjective = 0;
	}
	else {
		cout << "Error! solve_1_prototype_for_cluster_with_deleted_vertex\n";
	}
}

/*  
	solve 1-opt for cluster with inserted vertex
	inserted vertex is currentPrototype.assignedVertices.items[currentPrototype.assignedVertices.size-1]
*/
void solve_1_prototype_for_cluster_with_inserted_vertex(Prototype& currentPrototype) {
	int insertedVertexIndex = currentPrototype.assignedVertices.items[currentPrototype.assignedVertices.size - 1]; 
	double dist = 0;
	if ( currentPrototype.interior.edgeIndex != -1) {
		dist = min(shortestPathMatrix[insertedVertexIndex][edgeArray[currentPrototype.interior.edgeIndex].firstEndpointIndex] + currentPrototype.interior.length, shortestPathMatrix[insertedVertexIndex][edgeArray[currentPrototype.interior.edgeIndex].secondEndpointIndex] + edgeArray[currentPrototype.interior.edgeIndex].length - currentPrototype.interior.length);
	} else {
		dist = shortestPathMatrix[insertedVertexIndex][currentPrototype.vertexIndex];
	}
	BestOneClusterSolution bestSolution;
	bestSolution.objective = currentPrototype.clusterObjective + dist * dist;
	currentPrototype.clusterObjective = bestSolution.objective;
	int clusterSize = currentPrototype.assignedVertices.size;
	for (int i = 0; i < clusterSize; i++) {
		for (int j = 0; j < clusterSize; j++) {
			int currentEdge = isEdgeExist[currentPrototype.assignedVertices.items[i]][currentPrototype.assignedVertices.items[j]];
			if (currentEdge != -1) {
				double bottleneckPoint = 0.5*(edgeArray[currentEdge].length + shortestPathMatrix[insertedVertexIndex][edgeArray[currentEdge].secondEndpointIndex] - shortestPathMatrix[insertedVertexIndex][edgeArray[currentEdge].firstEndpointIndex]);
				double epsilon = shortestPathMatrix[insertedVertexIndex][edgeArray[currentEdge].firstEndpointIndex] + abs(bottleneckPoint) - dist;
				double intervalPoint1 = -1, intervalPoint2 = -1, intervalPoint3 = -1, intervalPoint4 = -1; 
				if (bottleneckPoint > 0 && bottleneckPoint < edgeArray[currentEdge].length) {
					intervalPoint1 = 0;
					intervalPoint2 = bottleneckPoint - epsilon;
					intervalPoint3 = bottleneckPoint + epsilon;
					intervalPoint4 = edgeArray[currentEdge].length;
				}
				else if (bottleneckPoint <= 0) {
					intervalPoint1 = max(0.0, bottleneckPoint + epsilon);
					intervalPoint2 = edgeArray[currentEdge].length;
				}
				else {
					intervalPoint1 = 0.0;
					intervalPoint2 = max(edgeArray[currentEdge].length, bottleneckPoint - epsilon);
				}
				solve_1_prototype_for_edge(currentEdge, currentPrototype, bestSolution, intervalPoint1, intervalPoint2);
				solve_1_prototype_for_edge(currentEdge, currentPrototype, bestSolution, intervalPoint3, intervalPoint4);
			}
		}
	}
	if (bestSolution.objective < currentPrototype.clusterObjective) {
		changePrototypeToBestOneClusterSolution(currentPrototype, bestSolution);
    }
}

void copyPrototypesAndUpdateTabu(Prototype* from, Prototype* to, bool* remainingPrototypes) {
	for (int cluster = 1; cluster <= clustersNumber; cluster++) {
		if (!remainingPrototypes[cluster]) {
			to[cluster].vertexIndex = from[cluster].vertexIndex;
			to[cluster].interior.edgeIndex = from[cluster].interior.edgeIndex;
			to[cluster].interior.length = from[cluster].interior.length;
            updateTabuInformation(to[cluster]);
		}
	}
}

void copyAllInformation(Clustering& from, Clustering& to) {
	for (int cluster = 1; cluster <= clustersNumber; cluster++) {
        to.prototypeArray[cluster].vertexIndex = from.prototypeArray[cluster].vertexIndex;
		to.prototypeArray[cluster].interior.edgeIndex = from.prototypeArray[cluster].interior.edgeIndex;
		to.prototypeArray[cluster].interior.length = from.prototypeArray[cluster].interior.length;
	}
    to.objective = from.objective;
}

void k_means(Clustering& currentSolution, bool* remainingPrototypes, char remove_degeneracyGen) {
	bool isDifferent = true;
	counter = 0;
	do {
		isDifferent = allocate(currentSolution);	
		if (isDifferent || counter == 0) {
			if (remove_degeneracy(remainingPrototypes, currentSolution.prototypeArray, remove_degeneracyGen)) {
				isDifferent = true;
                continue;
            }
			solve_1_prototype_for_every_cluster(currentSolution);
			isDifferent = true;
		}
		counter++;
        time (&endTime);
	} while(isDifferent && counter < 100 && difftime (endTime,startTime) < timeLimit);
    if (isDifferent) {
        for (int currentVertex = 1; currentVertex <= verticesNumber; currentVertex++) {
            vertexArray[currentVertex].distanceBetweenVertexAndPrototype = computeDistance(currentSolution.prototypeArray[vertexArray[currentVertex].clusterIndex], currentVertex);
        }
    }
}

void copyPrototypeWithoutAssignedVertices(Prototype& from, Prototype& to) {
	to.clusterObjective = from.clusterObjective;
	to.interior = from.interior;
	to.vertexIndex = from.vertexIndex;
}

void copyPrototypeWithAssignedVerticesSize(Prototype& from, Prototype& to) {
	copyPrototypeWithoutAssignedVertices(from, to);
	to.assignedVertices.size = from.assignedVertices.size;
}

void recomputeDistanceForCluster(Prototype& currentPrototype) {
	int currentVertex;
	for (int index = 0; index < currentPrototype.assignedVertices.size; index++) {
		currentVertex = currentPrototype.assignedVertices.items[index];
		vertexArray[currentVertex].distanceBetweenVertexAndPrototype = computeDistance(currentPrototype, currentVertex);
    }
}

void h_means(Clustering& currentSolution) {
    Prototype* prototypeArray = currentSolution.prototypeArray;
    Prototype oldPrototype1, oldPrototype2;
	int cluster1, cluster2;
	bool isFound;
	counter = 0;
	do {
		counter++;
		isFound = false;
        for (cluster1 = 1; cluster1 <= clustersNumber && !isFound; cluster1++) {
			int cluster1Size = prototypeArray[cluster1].assignedVertices.size;
			if (cluster1Size == 1) { continue; }
			copyPrototypeWithAssignedVerticesSize(prototypeArray[cluster1], oldPrototype1);
			for (int interchangeVertexIndex = 0; interchangeVertexIndex < cluster1Size; interchangeVertexIndex++) {
				int currentVertex = currentSolution.prototypeArray[cluster1].assignedVertices.items[interchangeVertexIndex];
				prototypeArray[cluster1].assignedVertices.items[interchangeVertexIndex] = prototypeArray[cluster1].assignedVertices.items[cluster1Size - 1];
				prototypeArray[cluster1].assignedVertices.size--;
				solve_1_prototype_for_cluster_with_deleted_vertex(prototypeArray[cluster1], currentVertex);
				for (cluster2 = 1; cluster2 <= clustersNumber; cluster2++) {
					if (cluster2 == cluster1) { continue; }
					int cluster2Size = prototypeArray[cluster2].assignedVertices.size;
					copyPrototypeWithAssignedVerticesSize(prototypeArray[cluster2], oldPrototype2);
					int vertexIndex = 0;
					for (; vertexIndex < cluster2Size; vertexIndex++) {
						if (isEdgeExist[prototypeArray[cluster2].assignedVertices.items[vertexIndex]][currentVertex] != -1) {
							break;
						}
					}
					if (vertexIndex == cluster2Size) {
						continue;
					}
					prototypeArray[cluster2].assignedVertices.items[prototypeArray[cluster2].assignedVertices.size++] = currentVertex;
					solve_1_prototype_for_cluster_with_inserted_vertex(prototypeArray[cluster2]);
					if (oldPrototype1.clusterObjective + oldPrototype2.clusterObjective > prototypeArray[cluster1].clusterObjective + prototypeArray[cluster2].clusterObjective) {
						currentSolution.objective -= oldPrototype1.clusterObjective + oldPrototype2.clusterObjective;
						currentSolution.objective += prototypeArray[cluster1].clusterObjective + prototypeArray[cluster2].clusterObjective;
						clearTabuInformation(oldPrototype1);
						clearTabuInformation(oldPrototype2);
						updateTabuInformation(prototypeArray[cluster1]);		
						updateTabuInformation(prototypeArray[cluster2]);
						vertexArray[currentVertex].clusterIndex = cluster2;
						recomputeDistanceForCluster(prototypeArray[cluster1]);
						recomputeDistanceForCluster(prototypeArray[cluster2]);
						isFound = true;
						break;
					} 
					copyPrototypeWithoutAssignedVertices(oldPrototype2, prototypeArray[cluster2]);
					prototypeArray[cluster2].assignedVertices.size--;
				}
				if (isFound) {
					break;
				}
				prototypeArray[cluster1].assignedVertices.items[cluster1Size - 1] = currentSolution.prototypeArray[cluster1].assignedVertices.items[interchangeVertexIndex];
				prototypeArray[cluster1].assignedVertices.items[interchangeVertexIndex] = currentVertex;
				prototypeArray[cluster1].assignedVertices.size++;
				copyPrototypeWithoutAssignedVertices(oldPrototype1, prototypeArray[cluster1]);
			}
		}
	} while (isFound && counter < 100);
}

void computeSecondMinDistanceAndCluster(Prototype* prototypeArray, int currentCluster, int currentVertex) {
    double minDistance = MAX_DOUBLE;
    int minCluster = -1;
	double curDistance;
    int cluster = 1;
    if (currentCluster > clustersNumber) {currentCluster = clustersNumber + 1;}
    for (; cluster < currentCluster; cluster++) {
        curDistance = computeDistance(prototypeArray[cluster], currentVertex); 
        if (curDistance < minDistance) {
            minDistance = curDistance;
			minCluster = cluster;
		}
    }
    for (cluster = currentCluster + 1; cluster <= clustersNumber; cluster++) {
        curDistance = computeDistance(prototypeArray[cluster], currentVertex); 
        if (curDistance < minDistance) {
            minDistance = curDistance;
			minCluster = cluster;
		}
    }
    vertexArray[currentVertex].secondMinimumDistance = minDistance;
	vertexArray[currentVertex].clusterIndexOfSecondMinimumDistance = minCluster;
}

double computeObjectiveWithNewPrototype(Prototype* prototypeArray, int changedCluster, Prototype& newPrototype) {
	double objective = 0;
	int currentVertex, clusterSize;
	double newDistanceBetweenProtAndVert;
	for (int currentCluster = 1; currentCluster <= clustersNumber; currentCluster++) {
		if (currentCluster == changedCluster) { continue; }
		clusterSize = prototypeArray[currentCluster].assignedVertices.size;
		for (int currentVertexIndex = 0; currentVertexIndex < clusterSize; currentVertexIndex++) {
			currentVertex = prototypeArray[currentCluster].assignedVertices.items[currentVertexIndex]; 
			newDistanceBetweenProtAndVert = computeDistance(newPrototype, currentVertex);
			if (vertexArray[currentVertex].distanceBetweenVertexAndPrototype < newDistanceBetweenProtAndVert) {
				objective += vertexArray[currentVertex].distanceBetweenVertexAndPrototype * vertexArray[currentVertex].distanceBetweenVertexAndPrototype;
			}
			else {
				objective += newDistanceBetweenProtAndVert * newDistanceBetweenProtAndVert;
			}
		}	
	}
	if (changedCluster > clustersNumber || changedCluster < 0) { return objective; }
	clusterSize = prototypeArray[changedCluster].assignedVertices.size;
	for (int currentVertexIndex = 0; currentVertexIndex < clusterSize; currentVertexIndex++) {
		currentVertex = prototypeArray[changedCluster].assignedVertices.items[currentVertexIndex];
		newDistanceBetweenProtAndVert = computeDistance(newPrototype, currentVertex);
		if (vertexArray[currentVertex].secondMinimumDistance < newDistanceBetweenProtAndVert) {
			objective += vertexArray[currentVertex].secondMinimumDistance * vertexArray[currentVertex].secondMinimumDistance;
		}
		else {
			objective += newDistanceBetweenProtAndVert * newDistanceBetweenProtAndVert;
		}
	}
	return objective;
}

double computeObjectiveWithNewVertexPrototype(Prototype* prototypeArray, int newPrototypeVertex, int changedCluster) {
	Prototype newPrototype;
	newPrototype.vertexIndex = newPrototypeVertex;
	newPrototype.interior.edgeIndex = -1;
	return computeObjectiveWithNewPrototype(prototypeArray, changedCluster, newPrototype);
}

double computeObjectiveWithNewInnerEdgePrototype(Prototype* prototypeArray, int newPrototypeEdge, double e, int changedCluster) {
	Prototype newPrototype;
	newPrototype.vertexIndex = -1;
	newPrototype.interior.edgeIndex = newPrototypeEdge;
	newPrototype.interior.length = e;
	return computeObjectiveWithNewPrototype(prototypeArray, changedCluster, newPrototype);
}

void recomputeObjectiveAndUpdateInfo(Clustering& currentSolution, Prototype& newPrototype, int changedCluster) {
	Prototype* prototypeArray = currentSolution.prototypeArray;
	double distance;
	currentSolution.objective = 0;
	for (int currentCluster = 1; currentCluster <= clustersNumber; currentCluster++) {
		prototypeArray[currentCluster].assignedVertices.size = 0;
		prototypeArray[currentCluster].clusterObjective = 0;
		if (currentCluster == changedCluster) {
			clearTabuInformation(prototypeArray[currentCluster]);
			prototypeArray[currentCluster].vertexIndex = newPrototype.vertexIndex;
			prototypeArray[currentCluster].interior.edgeIndex = newPrototype.interior.edgeIndex;
			prototypeArray[currentCluster].interior.length = newPrototype.interior.length;
			updateTabuInformation(prototypeArray[currentCluster]);
		}
	}
	int vertexCluster;
	for (int currentVertex = 1; currentVertex <= verticesNumber; currentVertex++) {
		if (vertexArray[currentVertex].distanceBetweenVertexAndPrototype > vertexArray[currentVertex].secondMinimumDistance) {
			vertexArray[currentVertex].clusterIndex = vertexArray[currentVertex].clusterIndexOfSecondMinimumDistance;
			vertexArray[currentVertex].distanceBetweenVertexAndPrototype = vertexArray[currentVertex].secondMinimumDistance;
			computeSecondMinDistanceAndCluster(prototypeArray, vertexArray[currentVertex].clusterIndex, currentVertex); 
		}
		distance = computeDistance(newPrototype, currentVertex);
		if (vertexArray[currentVertex].clusterIndex != changedCluster) {
			if (vertexArray[currentVertex].distanceBetweenVertexAndPrototype <= distance) {
				if (distance < vertexArray[currentVertex].secondMinimumDistance) {
					vertexArray[currentVertex].secondMinimumDistance = distance;
					vertexArray[currentVertex].clusterIndexOfSecondMinimumDistance = changedCluster;
				}
			}
			else {
				if (vertexArray[currentVertex].distanceBetweenVertexAndPrototype < vertexArray[currentVertex].secondMinimumDistance) {
					vertexArray[currentVertex].secondMinimumDistance = vertexArray[currentVertex].distanceBetweenVertexAndPrototype;
					vertexArray[currentVertex].clusterIndexOfSecondMinimumDistance = vertexArray[currentVertex].clusterIndex;
				}
				vertexArray[currentVertex].clusterIndex = changedCluster;
				vertexArray[currentVertex].distanceBetweenVertexAndPrototype = distance;
			}
		}
		else {
			if (vertexArray[currentVertex].secondMinimumDistance <= distance) {
				vertexArray[currentVertex].clusterIndex = vertexArray[currentVertex].clusterIndexOfSecondMinimumDistance;
				vertexArray[currentVertex].distanceBetweenVertexAndPrototype = vertexArray[currentVertex].secondMinimumDistance;
				computeSecondMinDistanceAndCluster(prototypeArray, vertexArray[currentVertex].clusterIndex, currentVertex); 
			}
			else {
				vertexArray[currentVertex].clusterIndex = changedCluster;
				vertexArray[currentVertex].distanceBetweenVertexAndPrototype = distance;
			}
		}
		vertexCluster = vertexArray[currentVertex].clusterIndex;
		prototypeArray[vertexCluster].assignedVertices.items[prototypeArray[vertexCluster].assignedVertices.size++] = currentVertex;
		prototypeArray[vertexCluster].clusterObjective += vertexArray[currentVertex].distanceBetweenVertexAndPrototype * vertexArray[currentVertex].distanceBetweenVertexAndPrototype;
		if (vertexArray[currentVertex].clusterIndexOfSecondMinimumDistance == changedCluster) {
			computeSecondMinDistanceAndCluster(prototypeArray, vertexArray[currentVertex].clusterIndex, currentVertex);	
		}
	}
	for (int currentCluster = 1; currentCluster <= clustersNumber; currentCluster++) {
		currentSolution.objective += prototypeArray[currentCluster].clusterObjective;
	}	
}

/*return true if successful*/
bool j_means(Clustering& currentSolution) {
	Prototype* prototypeArray = currentSolution.prototypeArray;
	int currentBestChangedCluster;
	double tmpObjFunction;
	double minNewObjFunction = currentSolution.objective;
	int currentVertex, currentVertexIndex, currentCluster;
    int lastChangedPrototype = -1;
	Prototype newPrototype;
	for (currentCluster = 1; currentCluster <= clustersNumber; currentCluster++) {
		for (currentVertexIndex = 0; currentVertexIndex < prototypeArray[currentCluster].assignedVertices.size; currentVertexIndex++) {
			computeSecondMinDistanceAndCluster(prototypeArray, currentCluster, prototypeArray[currentCluster].assignedVertices.items[currentVertexIndex]); 
		}	
	}
	counter = 0;
	do {
		counter++;
		currentBestChangedCluster = -1;
		minNewObjFunction = currentSolution.objective;
		newPrototype.vertexIndex = -1;
		newPrototype.interior.edgeIndex = -1;
		for (currentVertex = 1; currentVertex <= verticesNumber; currentVertex++) {
			if (vertexArray[currentVertex].tabu) {
				currentCluster = vertexArray[currentVertex].clusterIndex;
				if (currentSolution.prototypeArray[currentCluster].vertexIndex == -1) {
					int currentEdge = currentSolution.prototypeArray[currentCluster].interior.edgeIndex;
				/*	if (edgeArray[currentEdge].firstEndpointIndex != currentVertex && 
						edgeArray[currentEdge].secondEndpointIndex != currentVertex) {
						cout << "j-means error!\n";
					}*/
					tmpObjFunction = computeObjectiveWithNewVertexPrototype(prototypeArray, currentVertex, currentCluster);			
					if (tmpObjFunction < minNewObjFunction) {
						minNewObjFunction = tmpObjFunction;
						currentBestChangedCluster = currentCluster;
						newPrototype.vertexIndex = currentVertex;
                        currentVertex = verticesNumber;
					}
				}
				continue;
			}
			for (currentCluster = 1; currentCluster <= clustersNumber; currentCluster++) {
				if (currentCluster == lastChangedPrototype) {continue;}
				tmpObjFunction = computeObjectiveWithNewVertexPrototype(prototypeArray, currentVertex, currentCluster);			
				if (tmpObjFunction < minNewObjFunction) {
					minNewObjFunction = tmpObjFunction;
					currentBestChangedCluster = currentCluster;
					newPrototype.vertexIndex = currentVertex;
                    
                    currentVertex = verticesNumber;
                    break;
				}
			}
		}
		if (newPrototype.vertexIndex != -1) {
			recomputeObjectiveAndUpdateInfo(currentSolution, newPrototype, currentBestChangedCluster);
            lastChangedPrototype = currentBestChangedCluster;
		}
        time (&endTime);
	} while (newPrototype.vertexIndex != -1 && counter < 20 && difftime (endTime,startTime) < timeLimit);
    if (counter > 1) {
        return true;
    }
    return false;
}

double computeObjectiveForPotentialAssignedVertices(Prototype& newPrototype, int changedCluster) {
	int currentVertex;
	double currentSum = 0;
	double currentMin, distanceBetweenNewProtAndVertex;
	
	for (int vertexIndex = 0; vertexIndex < potentialAssignedVert.size; vertexIndex++) {
		currentVertex = potentialAssignedVert.items[vertexIndex];
		if (vertexArray[currentVertex].clusterIndex != changedCluster) {
			currentMin = vertexArray[currentVertex].distanceBetweenVertexAndPrototype;
		}
		else {
			currentMin = vertexArray[currentVertex].secondMinimumDistance;
		}
		distanceBetweenNewProtAndVertex = computeDistance(newPrototype, currentVertex);
		if (distanceBetweenNewProtAndVertex < currentMin) {
			currentSum += distanceBetweenNewProtAndVertex * distanceBetweenNewProtAndVertex;
		}
		else {
			currentSum += currentMin * currentMin;
		}
	}
	return currentSum;
}

double computeObjectiveWithInnerEdgePrototype(Prototype* prototypeArray, int newPrototypeEdge, int changedCluster, Prototype& newPrototype, double bestMinObj) {
	int firstEndpointIndex = edgeArray[newPrototypeEdge].firstEndpointIndex;
    int secondEndpointIndex = edgeArray[newPrototypeEdge].secondEndpointIndex;
	double intervalPoint1, intervalPoint2;
	if (vertexArray[firstEndpointIndex].clusterIndex == changedCluster) {
		intervalPoint2 = min(edgeArray[newPrototypeEdge].length, vertexArray[firstEndpointIndex].secondMinimumDistance);
	} else {
		intervalPoint2 = min(edgeArray[newPrototypeEdge].length, vertexArray[firstEndpointIndex].distanceBetweenVertexAndPrototype);
	}
	vertexArray[firstEndpointIndex].bottleneckPoint = 0.5*(edgeArray[newPrototypeEdge].length + shortestPathMatrix[firstEndpointIndex][secondEndpointIndex]);
	if (vertexArray[firstEndpointIndex].bottleneckPoint < intervalPoint2) {
		intervalPoint2 = vertexArray[firstEndpointIndex].bottleneckPoint;
	}		
	if (vertexArray[secondEndpointIndex].clusterIndex == changedCluster) {
		intervalPoint1 = max(0.0, edgeArray[newPrototypeEdge].length - vertexArray[secondEndpointIndex].secondMinimumDistance);
	} else {
		intervalPoint1 = max(0.0, edgeArray[newPrototypeEdge].length - vertexArray[secondEndpointIndex].distanceBetweenVertexAndPrototype);
	}
	vertexArray[secondEndpointIndex].bottleneckPoint = 0.5*(edgeArray[newPrototypeEdge].length - shortestPathMatrix[secondEndpointIndex][firstEndpointIndex]);
	if (vertexArray[secondEndpointIndex].bottleneckPoint > intervalPoint1) {
		intervalPoint1 = vertexArray[secondEndpointIndex].bottleneckPoint;
	}
	if (intervalPoint1 > intervalPoint2 || intervalPoint1 > edgeArray[newPrototypeEdge].length || intervalPoint2 < 0.0) {
        return MAX_DOUBLE;
    }
	if (intervalPoint1 == intervalPoint2) {
		newPrototype.vertexIndex = -1;
		newPrototype.interior.edgeIndex = newPrototypeEdge;
		newPrototype.interior.length = intervalPoint1;
		return computeObjectiveWithNewInnerEdgePrototype(prototypeArray, newPrototypeEdge, intervalPoint1, changedCluster);
	}
	intervalArray.size = 0;
	double leftShift, rightShift;
	double currentMin;
	double sumOfNonAssignedVertices = 0;
	leftShift = intervalPoint1;
	rightShift = edgeArray[newPrototypeEdge].length - intervalPoint2;
	shPath.sum_pathThroughFirstEndpoint = 0;
	shPath.sum_pathThroughSecondEndpoint = 0;
	shPath.pathThroughFirstEndpoint.size = 0;
	shPath.pathThroughSecondEndpoint.size = 0;
	potentialAssignedVert.size = 0;    
    for (int currentVertex = 1; currentVertex <= verticesNumber; currentVertex++) {
		if (vertexArray[currentVertex].clusterIndex != changedCluster) {
			currentMin = vertexArray[currentVertex].distanceBetweenVertexAndPrototype;
		}
		else {
			currentMin = vertexArray[currentVertex].secondMinimumDistance;
		}
		if (shortestPathMatrix[currentVertex][firstEndpointIndex] + leftShift < currentMin || 
			shortestPathMatrix[currentVertex][secondEndpointIndex] + rightShift < currentMin) {
				vertexArray[currentVertex].bottleneckPoint = 0.5*(edgeArray[newPrototypeEdge].length + shortestPathMatrix[currentVertex][secondEndpointIndex] - shortestPathMatrix[currentVertex][firstEndpointIndex]);
				if (shortestPathMatrix[currentVertex][firstEndpointIndex] + leftShift < currentMin && vertexArray[currentVertex].bottleneckPoint > intervalPoint1) {
					intervalArray.items[intervalArray.size].endpoint = min(currentMin - shortestPathMatrix[currentVertex][firstEndpointIndex], intervalPoint2);  
					intervalArray.items[intervalArray.size].isLeftEndpoint = true;
					if (intervalPoint1 < vertexArray[currentVertex].bottleneckPoint && vertexArray[currentVertex].bottleneckPoint < intervalArray.items[intervalArray.size].endpoint) {
						intervalArray.items[intervalArray.size].endpoint = vertexArray[currentVertex].bottleneckPoint;
					}
					shPath.sum_pathThroughFirstEndpoint += shortestPathMatrix[currentVertex][firstEndpointIndex]; 	
					shPath.pathThroughFirstEndpoint.size++;
					intervalArray.items[intervalArray.size++].vertex = currentVertex;
				}
				if (shortestPathMatrix[currentVertex][secondEndpointIndex] + rightShift < currentMin && vertexArray[currentVertex].bottleneckPoint < intervalPoint2) {
					intervalArray.items[intervalArray.size].endpoint = max(intervalPoint1, intervalPoint2 - (currentMin - shortestPathMatrix[currentVertex][secondEndpointIndex]));  
					intervalArray.items[intervalArray.size].isLeftEndpoint = false;
					if (intervalArray.items[intervalArray.size].endpoint < vertexArray[currentVertex].bottleneckPoint && vertexArray[currentVertex].bottleneckPoint < intervalPoint2) {
						intervalArray.items[intervalArray.size].endpoint = vertexArray[currentVertex].bottleneckPoint;
					}
					intervalArray.items[intervalArray.size++].vertex = currentVertex;
				}
				potentialAssignedVert.items[potentialAssignedVert.size++] = currentVertex;
        }
		else {
			sumOfNonAssignedVertices += currentMin * currentMin;
			if (bestMinObj <= sumOfNonAssignedVertices) {
				break;
			}
		}
    }
	if (bestMinObj <= sumOfNonAssignedVertices) {
		return MAX_DOUBLE;
	}
	qsort(intervalArray.items, intervalArray.size, sizeof(Interval), compareInterval);
	
	double curFirstEndpoint = intervalPoint1, curSecondEndpoint;
	double leftSum, rightSum;
	Prototype tmpPrototype;
	tmpPrototype.vertexIndex = -1;
	tmpPrototype.interior.edgeIndex = newPrototypeEdge;
	newPrototype.vertexIndex = -1;
	newPrototype.interior.edgeIndex = newPrototypeEdge;
	double currentMinObjForAsVert = MAX_DOUBLE; 
	double tmpObjForAsVert = MAX_DOUBLE;
	int endpointIndex = 0;
    
    for (; endpointIndex < intervalArray.size; endpointIndex++) {
		curSecondEndpoint = intervalArray.items[endpointIndex].endpoint;
		if (curFirstEndpoint != curSecondEndpoint) {
			leftSum = shPath.sum_pathThroughFirstEndpoint + curFirstEndpoint * shPath.pathThroughFirstEndpoint.size;
			rightSum = shPath.sum_pathThroughSecondEndpoint + (edgeArray[newPrototypeEdge].length - curFirstEndpoint) * shPath.pathThroughSecondEndpoint.size;
			if (leftSum < rightSum) {
				if (curFirstEndpoint + (rightSum - leftSum) / (shPath.pathThroughFirstEndpoint.size + shPath.pathThroughSecondEndpoint.size) <= curSecondEndpoint) {
					tmpPrototype.interior.length = curFirstEndpoint + (rightSum - leftSum) / (shPath.pathThroughFirstEndpoint.size + shPath.pathThroughSecondEndpoint.size);
					tmpObjForAsVert = computeObjectiveForPotentialAssignedVertices(tmpPrototype, changedCluster);
				} else if (curSecondEndpoint == intervalPoint2 && curSecondEndpoint < edgeArray[newPrototypeEdge].length) {
                    tmpPrototype.interior.length = curSecondEndpoint;
                    tmpObjForAsVert = computeObjectiveForPotentialAssignedVertices(tmpPrototype, changedCluster);
                }
			}
            else if (curFirstEndpoint == intervalPoint1 && curFirstEndpoint > 0.0) {
                tmpPrototype.interior.length = curFirstEndpoint;
                tmpObjForAsVert = computeObjectiveForPotentialAssignedVertices(tmpPrototype, changedCluster);
            }
            
            if (tmpObjForAsVert < currentMinObjForAsVert) {
                currentMinObjForAsVert = tmpObjForAsVert;
                newPrototype.interior.length = tmpPrototype.interior.length;
            }
		}
		if (!intervalArray.items[endpointIndex].isLeftEndpoint) {
			shPath.sum_pathThroughSecondEndpoint += shortestPathMatrix[intervalArray.items[endpointIndex].vertex][secondEndpointIndex]; 	
			shPath.pathThroughSecondEndpoint.size++;
		} else {
			shPath.sum_pathThroughFirstEndpoint -= shortestPathMatrix[intervalArray.items[endpointIndex].vertex][firstEndpointIndex]; 	
			shPath.pathThroughFirstEndpoint.size--;
		}
		curFirstEndpoint = curSecondEndpoint;
	}

	return sumOfNonAssignedVertices + currentMinObjForAsVert;
}

/*return true if successful*/
bool best_inner_prot_it(Clustering& currentSolution) {
	Prototype* prototypeArray = currentSolution.prototypeArray;
	int currentBestNewPrototypeVertex;
	int currentBestChangedCluster;
	int lastChangedPrototype = -1;
	double tmpObjFunction;
	double minNewObjFunction = currentSolution.objective;
	Prototype newPrototype;
	Prototype currentBestPrototype;
	int currentVertex, currentVertexIndex, currentCluster, currentEdge;
	for (currentCluster = 1; currentCluster <= clustersNumber; currentCluster++) {
		for (currentVertexIndex = 0; currentVertexIndex < prototypeArray[currentCluster].assignedVertices.size; currentVertexIndex++) {
			computeSecondMinDistanceAndCluster(prototypeArray, currentCluster, prototypeArray[currentCluster].assignedVertices.items[currentVertexIndex]); 
		}	
	}
	counter = 0;
	do {
		counter++;
		currentBestNewPrototypeVertex = -1;
		currentBestChangedCluster = -1;
		minNewObjFunction = currentSolution.objective;
		currentBestPrototype.vertexIndex = -1;
		currentBestPrototype.interior.edgeIndex = -1;
        for (currentVertex = 1; currentVertex <= verticesNumber; currentVertex++) {
			if (vertexArray[currentVertex].tabu) {
				currentCluster = vertexArray[currentVertex].clusterIndex;
				if (currentSolution.prototypeArray[currentCluster].vertexIndex == -1) {
					int currentEdge = currentSolution.prototypeArray[currentCluster].interior.edgeIndex;
					tmpObjFunction = computeObjectiveWithNewVertexPrototype(prototypeArray, currentVertex, currentCluster);			
					if (tmpObjFunction < minNewObjFunction) {
						minNewObjFunction = tmpObjFunction;
						currentBestChangedCluster = currentCluster;
						newPrototype.vertexIndex = currentVertex;
					}
				}
				continue;
			}
			for (currentCluster = 1; currentCluster <= clustersNumber; currentCluster++) {
				if (currentCluster == lastChangedPrototype) {continue;}
				tmpObjFunction = computeObjectiveWithNewVertexPrototype(prototypeArray, currentVertex, currentCluster);			
				if (tmpObjFunction < minNewObjFunction) {
					minNewObjFunction = tmpObjFunction;
					currentBestChangedCluster = currentCluster;
					newPrototype.vertexIndex = currentVertex;
				}
			}
		}
        for (currentEdge = 1; currentEdge <= edgesNumber; currentEdge++) {
			if (vertexArray[edgeArray[currentEdge].firstEndpointIndex].tabu || vertexArray[edgeArray[currentEdge].secondEndpointIndex].tabu) {
				currentCluster = -1;
				if (vertexArray[edgeArray[currentEdge].firstEndpointIndex].tabu && !vertexArray[edgeArray[currentEdge].secondEndpointIndex].tabu) {
					currentCluster = vertexArray[edgeArray[currentEdge].firstEndpointIndex].clusterIndex; 
				} else if (!vertexArray[edgeArray[currentEdge].firstEndpointIndex].tabu && vertexArray[edgeArray[currentEdge].secondEndpointIndex].tabu) {
					currentCluster = vertexArray[edgeArray[currentEdge].secondEndpointIndex].clusterIndex; 
				} else if (prototypeArray[currentCluster].interior.edgeIndex == currentEdge) {
					currentCluster = vertexArray[edgeArray[currentEdge].firstEndpointIndex].clusterIndex;
				}
				if (currentCluster != -1) {
					tmpObjFunction = computeObjectiveWithInnerEdgePrototype(prototypeArray, currentEdge, currentCluster, newPrototype, minNewObjFunction);
					if (tmpObjFunction < minNewObjFunction) {
						currentBestPrototype.vertexIndex = newPrototype.vertexIndex;
						currentBestPrototype.interior.edgeIndex = newPrototype.interior.edgeIndex;
						currentBestPrototype.interior.length = newPrototype.interior.length;
					
						minNewObjFunction = tmpObjFunction;
						currentBestChangedCluster = currentCluster;
					}
				}
				continue;
			}
			for (currentCluster = 1; currentCluster <= clustersNumber; currentCluster++) {
				if (currentCluster == lastChangedPrototype) {continue;}
				tmpObjFunction = computeObjectiveWithInnerEdgePrototype(prototypeArray, currentEdge, currentCluster, newPrototype, minNewObjFunction);
				if (tmpObjFunction < minNewObjFunction) {
					currentBestPrototype.vertexIndex = newPrototype.vertexIndex;
					currentBestPrototype.interior.edgeIndex = newPrototype.interior.edgeIndex;
					currentBestPrototype.interior.length = newPrototype.interior.length;
					
					minNewObjFunction = tmpObjFunction;
					currentBestChangedCluster = currentCluster;
				}
			}
		}
        
		if (currentBestPrototype.interior.edgeIndex != -1 || currentBestPrototype.vertexIndex != -1) {
			recomputeObjectiveAndUpdateInfo(currentSolution, currentBestPrototype, currentBestChangedCluster);
			lastChangedPrototype = currentBestChangedCluster;
		} 
	} while ((currentBestPrototype.interior.edgeIndex != -1 || currentBestPrototype.vertexIndex != -1) && counter < 20);
	return true;
}

/*return true if successful*/
bool i_means(Clustering& currentSolution) {
	Prototype* prototypeArray = currentSolution.prototypeArray;
	int currentBestNewPrototypeVertex;
	int currentBestChangedCluster;
	int lastChangedPrototype = -1;
	double tmpObjFunction;
	double minNewObjFunction = currentSolution.objective;
	Prototype newPrototype;
	Prototype currentBestPrototype;
	int currentVertexIndex, currentCluster, currentEdge;
	for (currentCluster = 1; currentCluster <= clustersNumber; currentCluster++) {
		for (currentVertexIndex = 0; currentVertexIndex < prototypeArray[currentCluster].assignedVertices.size; currentVertexIndex++) {
			computeSecondMinDistanceAndCluster(prototypeArray, currentCluster, prototypeArray[currentCluster].assignedVertices.items[currentVertexIndex]); 
		}	
	}
	counter = 0;
	do {
		counter++;
		currentBestNewPrototypeVertex = -1;
		currentBestChangedCluster = -1;
		minNewObjFunction = currentSolution.objective;
		currentBestPrototype.vertexIndex = -1;
		currentBestPrototype.interior.edgeIndex = -1;
		for (currentEdge = 1; currentEdge <= edgesNumber; currentEdge++) {
			if (vertexArray[edgeArray[currentEdge].firstEndpointIndex].tabu || vertexArray[edgeArray[currentEdge].secondEndpointIndex].tabu) {
				currentCluster = -1;
				if (vertexArray[edgeArray[currentEdge].firstEndpointIndex].tabu && !vertexArray[edgeArray[currentEdge].secondEndpointIndex].tabu) {
					currentCluster = vertexArray[edgeArray[currentEdge].firstEndpointIndex].clusterIndex; 
				} else if (!vertexArray[edgeArray[currentEdge].firstEndpointIndex].tabu && vertexArray[edgeArray[currentEdge].secondEndpointIndex].tabu) {
					currentCluster = vertexArray[edgeArray[currentEdge].secondEndpointIndex].clusterIndex; 
				} else if (prototypeArray[currentCluster].interior.edgeIndex == currentEdge) {
					currentCluster = vertexArray[edgeArray[currentEdge].firstEndpointIndex].clusterIndex;
					if (vertexArray[edgeArray[currentEdge].firstEndpointIndex].clusterIndex != vertexArray[edgeArray[currentEdge].secondEndpointIndex].clusterIndex) {
						cout << "i-means error!\n";
					}
				}
				if (currentCluster != -1) {
					tmpObjFunction = computeObjectiveWithInnerEdgePrototype(prototypeArray, currentEdge, currentCluster, newPrototype, minNewObjFunction);
					if (tmpObjFunction < minNewObjFunction) {
						currentBestPrototype.vertexIndex = newPrototype.vertexIndex;
						currentBestPrototype.interior.edgeIndex = newPrototype.interior.edgeIndex;
						currentBestPrototype.interior.length = newPrototype.interior.length;
					
						minNewObjFunction = tmpObjFunction;
						currentBestChangedCluster = currentCluster;

                        currentEdge = edgesNumber;
					}
				}
				continue;
			}
			for (currentCluster = 1; currentCluster <= clustersNumber; currentCluster++) {
				if (currentCluster == lastChangedPrototype) {continue;}
				tmpObjFunction = computeObjectiveWithInnerEdgePrototype(prototypeArray, currentEdge, currentCluster, newPrototype, minNewObjFunction);
				if (tmpObjFunction < minNewObjFunction) {
					currentBestPrototype.vertexIndex = newPrototype.vertexIndex;
					currentBestPrototype.interior.edgeIndex = newPrototype.interior.edgeIndex;
					currentBestPrototype.interior.length = newPrototype.interior.length;
					
					minNewObjFunction = tmpObjFunction;
					currentBestChangedCluster = currentCluster;
                    currentEdge = edgesNumber;
                    break;
				}
			}
		}
		if (currentBestPrototype.interior.edgeIndex != -1 || currentBestPrototype.vertexIndex != -1) {
			recomputeObjectiveAndUpdateInfo(currentSolution, currentBestPrototype, currentBestChangedCluster);
			lastChangedPrototype = currentBestChangedCluster;
		}
        time (&endTime);
    } while ((currentBestPrototype.interior.edgeIndex != -1 || currentBestPrototype.vertexIndex != -1) && counter < 20 && difftime (endTime,startTime) < timeLimit);
    if (counter > 1) {
        return true;
    }
    return false;
}

void VNS_old(char * argv[]) {
    bestSolution.objective = MAX_DOUBLE;
    cout << "---------------\nStarting! VNS_old\n";
    timeLimit = atoi(argv[2]);
    time (&startTime); 
	bool* remainingPrototypes = new bool [clustersNumber + 1];
	for (int cluster = 1; cluster <= clustersNumber; cluster++) {
		remainingPrototypes[cluster] = false;
	}
	vertexDiscreteDistribution = new std::uniform_int_distribution<int> (1, verticesNumber);
	edgeDiscreteDistribution = new std::uniform_int_distribution<int> (1, edgesNumber);
	
    generate_initialSolution(argv[3][0], bestSolution.prototypeArray);
    k_means(bestSolution, remainingPrototypes, argv[5][0]);
	printPrototypes(bestSolution.prototypeArray);
	printf("Initial solution result after k-means: %10.4f\n", bestSolution.objective);

    int k = 1;
    time (&endTime);
	while ( difftime (endTime,startTime) < timeLimit) {
        copyPrototypesAndUpdateTabu(bestSolution.prototypeArray, currentSolution.prototypeArray, remainingPrototypes);
		shaking(remainingPrototypes, k, currentSolution.prototypeArray, argv[4][0]);
		k_means(currentSolution, remainingPrototypes, argv[5][0]);
		if (currentSolution.objective < bestSolution.objective) {
			copyAllInformation(currentSolution, bestSolution);	
			k = 1;
		}
		else {
            k++;
            if (k > clustersNumber ) {
                k = 1;
            }
            clearTabuInformationForNotRemaining(remainingPrototypes, currentSolution.prototypeArray);
        }
        time (&endTime);
    }
    clearAllTabuInformation();
    printPrototypes(bestSolution.prototypeArray);
	printf("Best Solution: %10.4f\n", bestSolution.objective);
    cout << "Time: " << difftime (endTime,startTime) << " sec\n";
	cout << "Finishing!\n---------------\n";
    copyAllInformation(bestSolution, currentSolution);	
    allocate(currentSolution);
    printf("Final allocate: %10.4f\n", currentSolution.objective);
}

void VND_2(char * argv[], Clustering& currentSolution, bool* remainingPrototypes) {
    bool k_means_call = false;
    int iteration = 0;
    do {
        k_means(currentSolution, remainingPrototypes, argv[5][0]);
        k_means_call = j_means(currentSolution);
        iteration++;
    } while (k_means_call && iteration < 20);
}

void VND_3(char * argv[], Clustering& currentSolution, bool* remainingPrototypes) {
    bool k_means_call = false;
    int iteration = 0;
    do {
        k_means(currentSolution, remainingPrototypes, argv[5][0]);
        k_means_call = j_means(currentSolution);
        iteration++;
    } while (k_means_call && iteration < 20);
    i_means(currentSolution);
}

void VND_k_ij(char * argv[], Clustering& currentSolution, bool* remainingPrototypes) {
    k_means(currentSolution, remainingPrototypes, argv[5][0]);
    best_inner_prot_it(currentSolution);    
}

void VNS_2(char * argv[]) {
    bestSolution.objective = MAX_DOUBLE;
    cout << "---------------\nStarting! VNS_2\n";
    timeLimit = atoi(argv[2]);
    time (&startTime); 
	bool* remainingPrototypes = new bool [clustersNumber + 1];
	for (int cluster = 1; cluster <= clustersNumber; cluster++) {
		remainingPrototypes[cluster] = false;
	}
	vertexDiscreteDistribution = new std::uniform_int_distribution<int> (1, verticesNumber);
	edgeDiscreteDistribution = new std::uniform_int_distribution<int> (1, edgesNumber);
	
    generate_initialSolution(argv[3][0], bestSolution.prototypeArray);
    VND_2(argv, bestSolution, remainingPrototypes); 
    printPrototypes(bestSolution.prototypeArray);
	printf("Initial solution result after VND_2: %10.4f\n", bestSolution.objective);

    int k = 1;
    time (&endTime);
	while ( difftime (endTime,startTime) < timeLimit) {
        copyPrototypesAndUpdateTabu(bestSolution.prototypeArray, currentSolution.prototypeArray, remainingPrototypes);
		shaking(remainingPrototypes, k, currentSolution.prototypeArray, argv[4][0]);
        VND_2(argv, currentSolution, remainingPrototypes); 
        if (currentSolution.objective < bestSolution.objective) {
			copyAllInformation(currentSolution, bestSolution);	
			k = 1;
		}
		else {
            k++;
            if (k > clustersNumber ) {
                k = 1;
            }
            clearTabuInformationForNotRemaining(remainingPrototypes, currentSolution.prototypeArray);
        }
        time (&endTime);
    }
    clearAllTabuInformation();
    printPrototypes(bestSolution.prototypeArray);
	printf("Best Solution: %10.4f\n", bestSolution.objective);
    cout << "Time: " << difftime (endTime,startTime) << " sec\n";
	cout << "Finishing!\n---------------\n";
    copyAllInformation(bestSolution, currentSolution);	
    allocate(currentSolution);
    printf("Final allocate: %10.4f\n", currentSolution.objective);
}

void printClusters(Clustering& solution)
{
   Prototype* prototypeArray = solution.prototypeArray;
   int* objectsClusters = new int[verticesNumber + 1];
   for (int cluster = 1; cluster <= clustersNumber; cluster++) {
      for (int vertexIndex = 0; vertexIndex < prototypeArray[cluster].assignedVertices.size; ++vertexIndex)
         objectsClusters[prototypeArray[cluster].assignedVertices.items[vertexIndex]] = cluster;
   }
   if (verticesNumber > 0)
   {
      cout << objectsClusters[1];
      for (int vertexIndex = 2; vertexIndex <= verticesNumber; ++vertexIndex)
         cout << ", " << objectsClusters[vertexIndex];
      cout << '\n';
   }
}

void VNS_3(char * argv[]) {
    bestSolution.objective = MAX_DOUBLE;
    cout << "---------------\nStarting! VNS_3\n";
    timeLimit = atoi(argv[2]);
    time (&startTime); 
	bool* remainingPrototypes = new bool [clustersNumber + 1];
	for (int cluster = 1; cluster <= clustersNumber; cluster++) {
		remainingPrototypes[cluster] = false;
	}
	vertexDiscreteDistribution = new std::uniform_int_distribution<int> (1, verticesNumber);
	edgeDiscreteDistribution = new std::uniform_int_distribution<int> (1, edgesNumber);
	
    generate_initialSolution(argv[3][0], bestSolution.prototypeArray);
    VND_3(argv, bestSolution, remainingPrototypes); 
    printPrototypes(bestSolution.prototypeArray);
    printf("Initial solution result after VND_3: %10.4f\n", bestSolution.objective);

    int k = 1;
    time (&endTime);
	while ( difftime (endTime,startTime) < timeLimit) {
        copyPrototypesAndUpdateTabu(bestSolution.prototypeArray, currentSolution.prototypeArray, remainingPrototypes);
		shaking(remainingPrototypes, k, currentSolution.prototypeArray, argv[4][0]);
        VND_3(argv, currentSolution, remainingPrototypes); 
        if (currentSolution.objective < bestSolution.objective) {
			copyAllInformation(currentSolution, bestSolution);	
			k = 1;
		}
		else {
            k++;
            if (k > clustersNumber ) {
                k = 1;
            }
            clearTabuInformationForNotRemaining(remainingPrototypes, currentSolution.prototypeArray);
        }
        time (&endTime);
    }
    clearAllTabuInformation();
    printPrototypes(bestSolution.prototypeArray);
	printf("Best Solution: %10.4f\n", bestSolution.objective);
   printClusters(bestSolution);
    cout << "Time: " << difftime (endTime,startTime) << " sec\n";
	cout << "Finishing!\n---------------\n";
    copyAllInformation(bestSolution, currentSolution);	
    allocate(currentSolution);
    printf("Final allocate: %10.4f\n", currentSolution.objective);
}

void VNS_kji(char * argv[]) {
    bestSolution.objective = MAX_DOUBLE;
    cout << "---------------\nStarting! VNS_kji\n";
    timeLimit = atoi(argv[2]);
    time (&startTime); 
	bool* remainingPrototypes = new bool [clustersNumber + 1];
	for (int cluster = 1; cluster <= clustersNumber; cluster++) {
		remainingPrototypes[cluster] = false;
	}
	vertexDiscreteDistribution = new std::uniform_int_distribution<int> (1, verticesNumber);
	edgeDiscreteDistribution = new std::uniform_int_distribution<int> (1, edgesNumber);
	
    generate_initialSolution(argv[3][0], bestSolution.prototypeArray);
    VND_k_ij(argv, bestSolution, remainingPrototypes);
    printPrototypes(bestSolution.prototypeArray);
    printf("Initial solution result after VND_kji: %10.4f\n", bestSolution.objective);

    int k = 1;
    time (&endTime);
	while ( difftime (endTime,startTime) < timeLimit) {
        copyPrototypesAndUpdateTabu(bestSolution.prototypeArray, currentSolution.prototypeArray, remainingPrototypes);
		shaking(remainingPrototypes, k, currentSolution.prototypeArray, argv[4][0]);
        VND_k_ij(argv, currentSolution, remainingPrototypes);
        if (currentSolution.objective < bestSolution.objective) {
			copyAllInformation(currentSolution, bestSolution);	
			k = 1;
		}
		else {
            k++;
            if (k > clustersNumber ) {
                k = 1;
            }
            clearTabuInformationForNotRemaining(remainingPrototypes, currentSolution.prototypeArray);
        }
        time (&endTime);
    }
    clearAllTabuInformation();
    printPrototypes(bestSolution.prototypeArray);
	printf("Best Solution: %10.4f\n", bestSolution.objective);
    cout << "Time: " << difftime (endTime,startTime) << " sec\n";
	cout << "Finishing!\n---------------\n";
    copyAllInformation(bestSolution, currentSolution);	
    allocate(currentSolution);
    printf("Final allocate: %10.4f\n", currentSolution.objective);
}

void local_search_h_means(char * argv[]) {
    int timeLimit = atoi(argv[2]);
    time (&startTime); 
	bool* remainingPrototypes = new bool [clustersNumber + 1];
	for (int cluster = 1; cluster <= clustersNumber; cluster++) {
		remainingPrototypes[cluster] = false;
	}
	vertexDiscreteDistribution = new std::uniform_int_distribution<int> (1, verticesNumber);
	edgeDiscreteDistribution = new std::uniform_int_distribution<int> (1, edgesNumber);
	
    generate_initialSolution(argv[3][0], bestSolution.prototypeArray);
	allocate(bestSolution);
    h_means(bestSolution);
    clearTabuInformationForNotRemaining(remainingPrototypes, bestSolution.prototypeArray);
    
    long numberOfCalls = 1;
    long totalNumberOfSteps = counter;
    time (&endTime);
	while ( difftime (endTime,startTime) < timeLimit) {
        generate_initialSolution(argv[3][0], currentSolution.prototypeArray);
        allocate(currentSolution);
        h_means(currentSolution);
		if (currentSolution.objective < bestSolution.objective) {
			copyAllInformation(currentSolution, bestSolution);	
		}
		numberOfCalls++;
        totalNumberOfSteps += counter;
        clearTabuInformationForNotRemaining(remainingPrototypes, currentSolution.prototypeArray);
        time (&endTime);
    }
    printPrototypes(bestSolution.prototypeArray);
	printf("Best Solution: %10.4f\n", bestSolution.objective);
    cout << "Number of iterations: " << numberOfCalls << '\n';
    cout << "Total number of steps: " << totalNumberOfSteps << '\n';
    cout << "Time: " << difftime (endTime,startTime) << " sec\n";
	cout << "Finishing!\n---------------\n";
}

void local_search_k_means(char * argv[]) {
    int timeLimit = atoi(argv[2]);
    time (&startTime); 
	bool* remainingPrototypes = new bool [clustersNumber + 1];
	for (int cluster = 1; cluster <= clustersNumber; cluster++) {
		remainingPrototypes[cluster] = false;
	}
	vertexDiscreteDistribution = new std::uniform_int_distribution<int> (1, verticesNumber);
	edgeDiscreteDistribution = new std::uniform_int_distribution<int> (1, edgesNumber);
	
    generate_initialSolution(argv[3][0], bestSolution.prototypeArray);
	k_means(bestSolution, remainingPrototypes, argv[5][0]);
    clearTabuInformationForNotRemaining(remainingPrototypes, bestSolution.prototypeArray);
    
    long numberOfCalls = 1;
    long totalNumberOfSteps = counter;
    time (&endTime);
	while ( difftime (endTime,startTime) < timeLimit) {
        generate_initialSolution(argv[3][0], currentSolution.prototypeArray);
	    k_means(currentSolution, remainingPrototypes, argv[5][0]);
		if (currentSolution.objective < bestSolution.objective) {
			copyAllInformation(currentSolution, bestSolution);	
		}
		numberOfCalls++;
        totalNumberOfSteps += counter;
        clearTabuInformationForNotRemaining(remainingPrototypes, currentSolution.prototypeArray);
        time (&endTime);
    }
    printPrototypes(bestSolution.prototypeArray);
	printf("Best Solution: %10.4f\n", bestSolution.objective);
    cout << "Number of iterations: " << numberOfCalls << '\n';
    cout << "Total number of steps: " << totalNumberOfSteps << '\n';
    cout << "Time: " << difftime (endTime,startTime) << " sec\n";
	cout << "Finishing!\n---------------\n";
}

void local_search_j_means(char * argv[]) {
    int timeLimit = atoi(argv[2]);
    time (&startTime); 
	bool* remainingPrototypes = new bool [clustersNumber + 1];
	for (int cluster = 1; cluster <= clustersNumber; cluster++) {
		remainingPrototypes[cluster] = false;
	}
	vertexDiscreteDistribution = new std::uniform_int_distribution<int> (1, verticesNumber);
	edgeDiscreteDistribution = new std::uniform_int_distribution<int> (1, edgesNumber);
	
    generate_initialSolution(argv[3][0], bestSolution.prototypeArray);
	allocate(bestSolution);
    j_means(bestSolution);
    clearTabuInformationForNotRemaining(remainingPrototypes, bestSolution.prototypeArray);
    
    long numberOfCalls = 1;
    long totalNumberOfSteps = counter;
    time (&endTime);
	while ( difftime (endTime,startTime) < timeLimit) {
        generate_initialSolution(argv[3][0], currentSolution.prototypeArray);
        allocate(currentSolution);
        j_means(currentSolution);
		if (currentSolution.objective < bestSolution.objective) {
			copyAllInformation(currentSolution, bestSolution);	
		}
		numberOfCalls++;
        totalNumberOfSteps += counter;
        clearTabuInformationForNotRemaining(remainingPrototypes, currentSolution.prototypeArray);
        time (&endTime);
    }
    printPrototypes(bestSolution.prototypeArray);
	printf("Best Solution: %10.4f\n", bestSolution.objective);
    cout << "Number of iterations: " << numberOfCalls << '\n';
    cout << "Total number of steps: " << totalNumberOfSteps << '\n';
    cout << "Time: " << difftime (endTime,startTime) << " sec\n";
	cout << "Finishing!\n---------------\n";
}

void local_search_best_inner(char * argv[]) {
    int timeLimit = atoi(argv[2]);
    time (&startTime); 
	bool* remainingPrototypes = new bool [clustersNumber + 1];
	for (int cluster = 1; cluster <= clustersNumber; cluster++) {
		remainingPrototypes[cluster] = false;
	}
	vertexDiscreteDistribution = new std::uniform_int_distribution<int> (1, verticesNumber);
	edgeDiscreteDistribution = new std::uniform_int_distribution<int> (1, edgesNumber);
	
    generate_initialSolution(argv[3][0], bestSolution.prototypeArray);
	allocate(bestSolution);
    best_inner_prot_it(bestSolution);
    clearTabuInformationForNotRemaining(remainingPrototypes, bestSolution.prototypeArray);
    
    long numberOfCalls = 1;
    long totalNumberOfSteps = counter;
    time (&endTime);
	while ( difftime (endTime,startTime) < timeLimit) {
        generate_initialSolution(argv[3][0], currentSolution.prototypeArray);
        allocate(currentSolution);
        best_inner_prot_it(currentSolution);
		if (currentSolution.objective < bestSolution.objective) {
			copyAllInformation(currentSolution, bestSolution);	
		}
		numberOfCalls++;
        totalNumberOfSteps += counter;
        clearTabuInformationForNotRemaining(remainingPrototypes, currentSolution.prototypeArray);
        time (&endTime);
    }
    printPrototypes(bestSolution.prototypeArray);
	printf("Best Solution: %10.4f\n", bestSolution.objective);
    cout << "Number of iterations: " << numberOfCalls << '\n';
    cout << "Total number of steps: " << totalNumberOfSteps << '\n';
    cout << "Time: " << difftime (endTime,startTime) << " sec\n";
	cout << "Finishing!\n---------------\n";
}


int main(int argc, char * argv[]) 
{
	if (argc < 6) 
	{
		exitOnError("Usage: program.x <filename>", "");
	}
	readInstanceAndInitializeVariables(argv[1]);
    Floyd–Warshall_algorithm();
//    VNS_old(argv);
//    VNS_2(argv);
    VNS_3(argv);
//    VNS_kji(argv);
/*	for (int it = 1; it <= 5; it++) {
        local_search_k_means(argv);
        local_search_h_means(argv);
        local_search_j_means(argv);
	    local_search_best_inner(argv);
    }*/
    return 0;
}

