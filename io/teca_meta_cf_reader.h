#ifndef teca_algorithm_h
#define teca_algorithm_h

#include "teca_algorithm.h"

TECA_SHARED_OBJECT_FORWARD_DECL(teca_cf_reader)

class teca_meta_cf_reader_internals;
using p_teca_meta_cf_reader_internals = std::shared_ptr<teca_meta_cf_reader_internals>;

// interface to teca pipeline architecture. all sources/readers
// filters, sinks/writers will implement this interface
class teca_meta_cf_reader : public teca_algorithm
{
public:
    TECA_ALGORITHM_STATIC_NEW(teca_meta_cf_reader)
    TECA_ALGORITHM_DELETE_COPY_ASSIGN(teca_meta_cf_reader)
    TECA_ALGORITHM_CLASS_NAME(teca_meta_cf_reader)
    ~teca_meta_cf_reader();

    // report/initialize to/from Boost program options
    // objects.
    TECA_GET_ALGORITHM_PROPERTIES_DESCRIPTION()
    TECA_SET_ALGORITHM_PROPERTIES()

    // adds a reader to the collection and at the same time specifies
    // how it will be used.
    void add_reader(const std::string &key,
        const p_teca_cf_reader &reader,
        int provides_time, int provides_geometry,
        const std::set<std::string> &variables);

    // adds a reader to the collection.
    void add_reader(const std::string &key,
        const p_teca_cf_reader &reader);

    // sets the reader that provides the time axis
    void set_time_reader(const std::string &key);

    // sets the reader that provides the mesh geometry
    void set_geometry_reader(const std::string &key);

    // adds to the list of variables that a reader will provide
    void add_variable_reader(const std::string &key,
        const std::string &variable);

    // sets the list of variable that a reader will provide.
    void set_variable_reader(const std::string &key,
        const std::set<std::string> &variable)


protected:
    teca_meta_cf_reader();

    teca_metadata get_output_metadata(unsigned int port,
        const std::vector<teca_metadata> &input_md) override;

    const_p_teca_dataset execute(unsigned int port,
        const std::vector<const_p_teca_dataset> &input_data,
        const teca_metadata &request) override;

private:
    teca_algorithm_internals *internals;
};

#endif
