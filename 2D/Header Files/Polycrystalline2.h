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

	// Add external polycrystalline shell and some method for its parsing from inclInCrysesNum
	vector<unique_ptr<ShellEdge2>*> _shellEdges;
	vector<unique_ptr<ShellNode2>*> _shellNodes;

	list<unique_ptr<Simplex2>*> _freeSimplexes;
	vector<unique_ptr<Edge2>*> _freeEdges;
	vector<unique_ptr<Node2>*> _freeNodes;

public:
	list<Crystallite2*> crystallites;

	void GenerateNodesEvenly(double* const polycrysSizeAxis, size_t* const minNodesNumAxis);
	unique_ptr<Edge2>* FindEdge(Node2& node0, Node2& node1);
	void GenerateSimplexesFromNodes(size_t minNodesNumAxis_0);
	void GenerateUniformMesh();
	void FitNodesToShellNodes();
	//void FitNodesToShellEdges(); // DEPRECATED
	void FitNodesToShellEdges();
	void FitMeshToShells();
	void ErasePtrsToNullptrFromVectors();
	void DeleteExternalNodes();
	void DeleteFarNodes();
	//void DivideExtendedSimplexes();
	//void DivideCrossingEdges();
	//void MinMaxEdges(double& min, double& max);
	Vector2 ShiftToFitMesh(const Node2& node);
	Vector2 ShiftToLavEdges(const Node2& node);
	//Vector2 ShiftToEquilateralSimplexes(const Node2& node);
	void DistributeNodesEvenly();
	//void MakeSimplexesEquilateral();
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