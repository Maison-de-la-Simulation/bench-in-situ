#pragma once

#include "HydroParams.hpp"
#include "HydroProblem.hpp"
#include "HydroSolver.hpp"
#include "HydroTypes.hpp"

#include <iostream>
#include <memory>
#include <string>

namespace hydro
{

class Engine
{
public:
    Engine() = default;
    Engine(const Engine& x) = default;
    Engine(Engine&& x) = default;
    virtual ~Engine() = default;
    Engine& operator=(const Engine& x) = default;
    Engine& operator=(Engine&& x) = default;

    virtual void print_configuration(std::ostream& os) const = 0;
    virtual void run() const = 0;
};

struct EngineNd : Engine
{
    EngineNd(const std::string& filename);

    void print_configuration(std::ostream& os) const final;
    void run() const final;

    std::shared_ptr<Params>  params;
    std::shared_ptr<Problem> problem;
    std::shared_ptr<Solver>       solver;
};

struct EngineFactory
{
    static std::shared_ptr<Engine> New(const std::string& filename);
};

}
