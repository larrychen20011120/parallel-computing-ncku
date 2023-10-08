#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <list>
#include "mpi.h"

#define MAX_POINT_COUNT 12000
#define MAIN_PROCESS 0
#define CODE 0
#define MAX_VALUE 999999999

using namespace std;

// define point position type
struct Point {
    /* data */
    int x, y;
    int originIndex;

    // overload operator
    bool operator<(Point other) const {
      if(x == other.x) {
         return (y < other.y) ? true : false;
      } else if (x < other.x) {
        return true;
      } else {
        return false;
      }
    }
    Point& operator=(const Point& other) {
        x = other.x;
        y = other.y;
        originIndex = other.originIndex;
        return *this;
    }
};

bool isPointUp(Point prev, Point curr, Point next);
void copyArr(Point*, Point*, int);

int main(int argc, char *argv[]) {

    // problem's valriables
    int n;
    string filename;
    Point localPoints[MAX_POINT_COUNT];
    // list implemented stack
    // bi-directed
    list<int> prevIndex;

    // mpi variables
    int currentId;
    int processCount;

    // start mpi
    MPI_Init(&argc, &argv);
    // retrieve id
    MPI_Comm_rank(MPI_COMM_WORLD, &currentId);
    // get the total process count
    MPI_Comm_size(MPI_COMM_WORLD, &processCount);

    // read data
    if (currentId == MAIN_PROCESS) {
        cin >> filename;
        ifstream ifs(filename);
        ifs >> n;
        for (int i = 0; i < n; ++i) {
            ifs >> localPoints[i].x >> localPoints[i].y;
            localPoints[i].originIndex = i;
        }
        ifs.close();
    }
    
    int localCounts[processCount];
    for (int i = 0; i < processCount; i++) {
        localCounts[i] = n / processCount;
    }
    localCounts[processCount-1] += (n % processCount);

    MPI_Bcast(&n, 1, MPI_INT, MAIN_PROCESS, MPI_COMM_WORLD);
    MPI_Bcast(localCounts, processCount, MPI_INT, MAIN_PROCESS, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    
    if (currentId == MAIN_PROCESS) {
        // init all local points
        // send
        for (int id = 1; id < processCount; id++) {
            int start = id * n / processCount;
            MPI_Send(localPoints + start, localCounts[id]*sizeof(Point), MPI_BYTE, id,
                     CODE, MPI_COMM_WORLD);
        }

    } else {
        // receive
        MPI_Recv(localPoints, localCounts[currentId]*sizeof(Point), MPI_BYTE, MAIN_PROCESS, 
                 CODE, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // local sort by std sort
    sort(localPoints, localPoints + localCounts[currentId]);

    
    // merge from local sort
    int recvDist = 2, sendDist = 1;
    int partner;
    while (recvDist <= processCount) {
        if (currentId % recvDist == 0) {
            // receiver
            partner = currentId + sendDist;
            if (partner < n) {
                // have sender to recive
                Point temp[MAX_POINT_COUNT];
                Point recv[MAX_POINT_COUNT];
                MPI_Recv(recv, localCounts[partner]*sizeof(Point), MPI_BYTE, partner,
                        CODE, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                copyArr(temp, localPoints, localCounts[currentId]);
                // for the end case comparison
                temp[localCounts[currentId]].x = MAX_VALUE;
                recv[localCounts[partner]].x = MAX_VALUE;

                // start merging
                for (int i = 0, j = 0; i < localCounts[currentId] || j < localCounts[partner];) {
                    if (temp[i] < recv[j]) {
                        localPoints[i+j] = temp[i];
                        i++;
                    } else {
                        localPoints[i+j] = recv[j];
                        j++;
                    }
                }
                
                // update tables
                for (int id = 0; id < processCount; id++) {
                    if (id % recvDist == 0) {
                        partner = id + sendDist;
                        if (partner < n) {
                            localCounts[id] += localCounts[partner];
                        }
                    }
                }
            }
        } else if (currentId % sendDist == 0) {
            // sender
            partner = currentId - sendDist;
            MPI_Send(localPoints, localCounts[currentId]*sizeof(Point), MPI_BYTE, partner,
                        CODE, MPI_COMM_WORLD);
        }
        
        recvDist *= 2; sendDist *= 2;

        // update all localCount to global value

        MPI_Barrier(MPI_COMM_WORLD);
    }

    // broadcast all points to all processes
    MPI_Bcast(localPoints, n*sizeof(Point), MPI_BYTE, MAIN_PROCESS, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    if (processCount >= 2) {
        // do parallel
        if (currentId == MAIN_PROCESS) {
            // search upper edge
            prevIndex.push_back(0);
            prevIndex.push_back(1);
            for (int i = 2; i < n; i++) {
                while (prevIndex.size() >= 2) {
                    int curr = prevIndex.back();
                    prevIndex.pop_back();

                    if (isPointUp(localPoints[prevIndex.back()], localPoints[curr], localPoints[i])) {
                        prevIndex.push_back(curr);
                        break;
                    }
                }
                prevIndex.push_back(i);
            }
            list<int>::iterator iter;
            for (iter = prevIndex.begin(); iter != prevIndex.end(); ++iter) {
                cout << localPoints[(*iter)].originIndex+1 << " ";
            }

            int numOfElement;
            MPI_Recv(&numOfElement, 1, MPI_INT, 1,
                    CODE, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            
            for (int i = 0; i < numOfElement; i++) {
                int ans;
                MPI_Recv(&ans, 1, MPI_INT, 1, CODE, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                cout << ans << " ";
            }

        }
        if (currentId == 1) {
            // search lower edge
            prevIndex.push_back(n-1);
            prevIndex.push_back(n-2);
            for (int i = n-3; i >= 0; --i) {
                while (prevIndex.size() >= 2) {
                    int curr = prevIndex.back();
                    prevIndex.pop_back();

                    if (isPointUp(localPoints[prevIndex.back()], localPoints[curr], localPoints[i])) {
                        prevIndex.push_back(curr);
                        break;
                    }
                }
                prevIndex.push_back(i);
            }
            prevIndex.pop_back();
            prevIndex.pop_front();
            
            int numOfElement = prevIndex.size();
            // send how many elements in stack
            MPI_Send(&numOfElement, 1, MPI_INT, MAIN_PROCESS,
                     CODE, MPI_COMM_WORLD);
            list<int>::iterator iter;
            for (iter = prevIndex.begin(); iter != prevIndex.end(); ++iter) {
                int ans = localPoints[(*iter)].originIndex+1;
                MPI_Send(&ans, 1, MPI_INT, MAIN_PROCESS,
                         CODE, MPI_COMM_WORLD);
            }
        }


    } else { // can't parallel
    
        // search upper edge
        prevIndex.push_back(0);
        prevIndex.push_back(1);
        for (int i = 2; i < n; i++) {
            while (prevIndex.size() >= 2) {
                int curr = prevIndex.back();
                prevIndex.pop_back();

                if (isPointUp(localPoints[prevIndex.back()], localPoints[curr], localPoints[i])) {
                    prevIndex.push_back(curr);
                    break;
                }
            }
            prevIndex.push_back(i);
        }

        list<int>::iterator iter;
        for (iter = prevIndex.begin(); iter != prevIndex.end(); ++iter) {
            cout << localPoints[(*iter)].originIndex+1 << " ";
        }

        prevIndex.clear();
        // search lower edge
        prevIndex.push_back(n-1);
        prevIndex.push_back(n-2);
        for (int i = n-3; i >= 0; --i) {
            while (prevIndex.size() >= 2) {
                int curr = prevIndex.back();
                prevIndex.pop_back();

                if (isPointUp(localPoints[prevIndex.back()], localPoints[curr], localPoints[i])) {
                    prevIndex.push_back(curr);
                    break;
                }
            }
            prevIndex.push_back(i);
        }
        prevIndex.pop_back();
        prevIndex.pop_front();

        for (iter = prevIndex.begin(); iter != prevIndex.end(); ++iter) {
            cout << localPoints[(*iter)].originIndex+1 << " ";
        }

    }



    MPI_Finalize();

    return 0;
}

// detect the intermidiate point is upper than the previous or not
bool inline isPointUp(Point prev, Point curr, Point next) {
    double slopeNext, slopeCurrent;
    if ((next.x - prev.x) != 0) {
        slopeNext = ((double)(next.y - prev.y)) / (next.x - prev.x);
    } else {
        slopeNext = MAX_VALUE;
    }
    if ((curr.x - prev.x) != 0) {
        slopeCurrent = ((double)(curr.y - prev.y)) / (curr.x - prev.x);
    } else {
        slopeCurrent = MAX_VALUE;
    }
    return (slopeCurrent > slopeNext) ? true : false;
}

void inline copyArr(Point* dest, Point* src, int count) {
    for (int i = 0; i < count; i++) {
        dest[i] = src[i];
    }
}