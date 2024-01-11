#include "Hydro.hpp"
#include "HydroEngine.hpp"
#include "Print.hpp"

#include <cstdlib>
#include <exception>
#include <memory>
#include <string>

#if !defined(NDEBUG)
#include <cfenv>
#endif // !defined(NDEBUG)

#if defined(Euler_ENABLE_PDI)
#include <pdi.h>
#endif // defined(Euler_ENABLE_PDI)

int main(int argc, char** argv)
{
#if !defined(NDEBUG)
#if defined(__GNUG__)
    feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW);
#endif // defined(__GNUG__)
#endif // !defined(NDEBUG)

#if defined(Euler_ENABLE_PDI)
    if (argc < 3)
    {
        std::cout << "usage: " << argv[0] << " <path to the ini file> <path to the yaml file>\n";
        return EXIT_FAILURE;
    }
#else
    if (argc < 2)
    {
        std::cout << "usage: " << argv[0] << " <path to the ini file>\n";
        return EXIT_FAILURE;
    }
#endif // defined(Euler_ENABLE_PDI)

    hydro::initialize(argc, argv);

#if defined(Euler_ENABLE_PDI)
    PC_tree_t conf = PC_parse_path(argv[2]);
    PDI_init(PC_get(conf, ".pdi"));
#endif // defined(Euler_ENABLE_PDI)
    try
    {
        auto engine = hydro::EngineFactory::New(argv[1]);

        std::ostringstream oss;
        hydro::print_configuration(oss);
        engine->print_configuration(oss);
        Print() << oss.str();

        engine->run();
    }
    catch (const std::exception& e)
    {
        std::cerr << e.what() << std::endl;
    }

#if defined(Euler_ENABLE_PDI)
    PDI_finalize();
    PC_tree_destroy(&conf);
#endif // defined(Euler_ENABLE_PDI)

    hydro::finalize();

    return EXIT_SUCCESS;
}
