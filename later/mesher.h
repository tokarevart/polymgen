#pragma once
#include "polytope.h"
#include "genparams.h"


namespace pmg {

template <typename Polytope>
class mesher
{
public:
    using polytope_type = Polytope;

    void run(polytope_type::real_type preferredLength, genparams<Polytope> genParams = genparams<Polytope>());

    mesher(const polytope_type* polytope)
        : m_polytope(polytope) {}


private:
    const polytope_type* m_polytope;
};

template <typename Polytope>
class mesher<std::vector<Polytope*>>
{
public:
    using polytope_type = Polytope;

    void run(polytope_type::real_type preferredLength, genparams<Polytope> genParams = genparams<Polytope>());

    mesher(const polytope_type* polytope)
        : m_polytope(polytope) {}


private:
    const polytope_type* m_polytope;
};

} // namespace pmg