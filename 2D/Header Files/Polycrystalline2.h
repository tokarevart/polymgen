#pragma once
#include <fstream>
#include <list>
#include <vector>
#include <memory>
#include "AllClassInclusions.h"
#include "AllClassDefinitions.h"


class Polycrystalline2
{
private:
	double _minShellNodesCoor[2];
	double _maxShellNodesCoor[2];
	double _l_min, _l_av, _l_max;

	std::vector<std::unique_ptr<ShellEdge2>> _shellEdges;
	std::vector<std::unique_ptr<ShellEdge2>> _shellNodes;

	std::list<Simplex2*> _freeSimplexes;
	std::vector<Edge2*> _freeEdges;
	std::vector<Node2*> _freeNodes;

public:
	std::list<Crystallite2*> crystallites;

	void Debug();
	void GenerateFreeNodesEvenly(double* const polycrysSizeAxis, size_t* const minNodesNumAxis);
	Edge2* FindFreeEdge(const Node2& node0, const Node2& node1);
	void GenerateFreeSimplexesFromFreeNodes(size_t minNodesNumAxis_0);
	void GenerateFreeUniformMesh(); // There may will not be significant performance improvement. But if you want you can parallize 1-st step.
	void FitFreeNodesToShellNodes();
	void FitFreeNodesToShellEdges();
	void FitFreeMeshToShells(); // There may will not be significant performance improvement. But if you want you can parallize 1-st step.
	//void DeleteExternalNodes();
	//void FillDataWithRemainingDependences(); // 1-st step can be well parallelized by nodes.
	//void DistributeNodesEvenly();
	//void GenerateFiniteElementMesh(); // Executes other functions.

	// Use when simplexes list already filled.
	const bool IsContaining(const Crystallite2& crys) const;
	const bool IsContaining(const Simplex2& simp) const;
	const bool IsContaining(const Edge2& edge) const;
	const bool IsContaining(const Node2& node) const;

	// Calculate based on position in space.
	const bool Contains(const Simplex2& simp) const;
	const bool Contains(const Edge2& edge) const;
	const bool Contains(const Node2& node) const;

	void InputData(std::ifstream& shell_nodes, std::ifstream& cryses_edges, std::ifstream& gener_params);
	void InputData(const std::vector<double>& shell_nodes, const std::vector<size_t>& cryses_edges, const std::vector<size_t>& gener_params);
	
	void OutputData(std::ofstream& nodesData, std::ofstream& feNodesData) const;
	void OutputData(std::vector<double>& nodesData, std::vector<size_t>& feNodesData) const;

	Polycrystalline2();
	Polycrystalline2(std::ifstream& shell_nodes, std::ifstream& cryses_edges, std::ifstream& gener_params);
	Polycrystalline2(const std::vector<double>& shell_nodes, const std::vector<size_t>& cryses_edges, const std::vector<size_t>& gener_params);
	~Polycrystalline2();
};