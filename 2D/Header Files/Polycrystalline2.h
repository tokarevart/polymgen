#pragma once
#include <fstream>
#include <list>
#include <vector>
#include <memory>
#include "ClassDefinitions.h"
#include "ClassInclusions.h"

using std::ifstream;
using std::ofstream;
using std::unique_ptr;
using std::list;
using std::vector;

class Polycrystalline2
{
private:
	double _minShellNodesCoor[2];
	double _maxShellNodesCoor[2];
	double _l_min, _l_av, _l_max;

	vector<unique_ptr<ShellEdge2>*> _shellEdges;
	vector<unique_ptr<ShellNode2>*> _shellNodes;

	list<unique_ptr<Simplex2>*> _freeSimplexes;
	vector<unique_ptr<Edge2>*> _freeEdges;
	vector<unique_ptr<Node2>*> _freeNodes;

public:
	list<Crystallite2*> crystallites;

	void Debug();
	void GenerateFreeNodesEvenly(double* const polycrysSizeAxis, size_t* const minNodesNumAxis);
	unique_ptr<Edge2>* FindFreeEdge(Node2& node0, Node2& node1);
	void GenerateFreeSimplexesFromFreeNodes(size_t minNodesNumAxis_0);
	void GenerateFreeUniformMesh(); // There may will not be significant performance improvement. But if you want you can parallize 1-st step.
	void FitFreeNodesToShellNodes();
	void FitFreeNodesToShellEdges();
	void FitFreeMeshToShells(); // There may will not be significant performance improvement. But if you want you can parallize 1-st step.
	//void DeleteExternalNodes();
	//void FillDataWithRemainingDependences(); // 1-st step can be well parallelized by nodes.
	//void DistributeNodesEvenly();
	//void GenerateFiniteElementMesh(); // Executes other functions.
	
	// Calculate based on position in space.
	const bool Contains(const Node2& node) const;

	void InputData(ifstream& shell_nodes, ifstream& cryses_edges, ifstream& gener_params);
	void InputData(const vector<double>& shell_nodes, const vector<size_t>& cryses_edges, const vector<size_t>& gener_params);
	
	void OutputData(ofstream& nodesData, ofstream& feNodesData) const;
	void OutputData(vector<double>& nodesData, vector<size_t>& feNodesData) const;

	Polycrystalline2();
	Polycrystalline2(ifstream& shell_nodes, ifstream& cryses_edges, ifstream& gener_params);
	Polycrystalline2(const vector<double>& shell_nodes, const vector<size_t>& cryses_edges, const vector<size_t>& gener_params);
	~Polycrystalline2();
};