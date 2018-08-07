#pragma once
#include <fstream>
#include <list>
#include <vector>
#include "AllClassInclusions.h"
#include "AllClassDefinitions.h"


class Polycrystalline2
{
private:
	double minShellNodesCoor[2];
	double maxShellNodesCoor[2];
	double l_min, l_av, l_max;

public:
	std::list<Crystallite2*> crystallites;

	std::list<Simplex2*> freeSimplexes;
	std::vector<Edge2*> freeEdges;
	std::vector<Node2*> freeNodes;

	void AddCrystallite();                          // May be useless if link from crystallite to polycrystalline will not be used.
	void AddCrystallite(Crystallite2* const& crys); // May be useless if link from crystallite to polycrystalline will not be used.
	
	size_t GenerateFreeNodesEvenly(double* const& polycrysSizeAxis, size_t* const& minNodesNumAxis);
	void GenerateFreeSimplexesFromFreeNodes();
	void GenerateFreeUniformMesh(); // There may will not be significant performance improvement. But if you want you can parallize 1-st step.
	//void FitFreeMeshToShells(); // There may will not be significant performance improvement. But if you want you can parallize 1-st step.
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