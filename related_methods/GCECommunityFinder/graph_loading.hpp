#ifndef GRAPH_LOADING_H
#define GRAPH_LOADING_H

#include "graph_representation.hpp"
#include "Range.hpp"

namespace graph_loading {

void loadBloomGraphMMAP(SimpleIntGraph &bg, const char *fileName);
void loadSimpleIntGraphFromFile(SimpleIntGraph &bg, const char *fileName);
struct readEdgeInvalidLineInDataException : public exception { }; // TODO: Proper error 
const char *readEdge(const char *cur, int &l, int &r) throw (readEdgeInvalidLineInDataException); 
} // namespace graph_loading

#endif