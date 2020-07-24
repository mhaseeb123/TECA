#include "teca_profiler.h"
#include "teca_memory_profiler.h"
#include "teca_common.h"
#include "teca_config.h"

#include <sys/time.h>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <stdint.h>
#include <strings.h>
#include <cstdlib>
#include <cstdlib>
#include <cstdio>

#include <map>
#include <list>
#include <vector>
#include <iomanip>
#include <limits>
#include <unordered_map>
#include <mutex>

#if defined(TECA_HAS_TIMEMORY)
// timemory include
#include "timemory/timemory.hpp"

using namespace tim::component;

// user bundle
struct teca_time_tag 
{};

// custom user bundle
using timing_bundle_t = user_bundle<0, teca_time_tag>;

// add components to time profiler
using measurement_t   = tim::component_tuple_t<wall_clock, timing_bundle_t>;
using bundle_t        = timing_bundle_t;

namespace measure
{
// ----------------------------------------------------------------------------
// tuple access helper functions
template <typename Tuple, typename F, std::size_t ...Indices>
void for_each_impl(Tuple&& tuple, F&& f, std::index_sequence<Indices...>) {
    using swallow = int[];
    (void)swallow{1,
        (f(std::get<Indices>(std::forward<Tuple>(tuple))), void(), int{})...
    };
}

// ----------------------------------------------------------------------------
// tuple access helper functions
template <typename Tuple, typename F>
void apply_tuple(Tuple&& tuple, F&& f) {
    constexpr std::size_t N = std::tuple_size<std::remove_reference_t<Tuple>>::value;
    for_each_impl(std::forward<Tuple>(tuple), std::forward<F>(f), 
                  std::make_index_sequence<N>{});
}

namespace time
{
// is global bundle configured
static bool is_configured = false;

// ----------------------------------------------------------------------------
// print static component names onto the stringstream
template<typename Tp>
void write_names(Tp comp, std::ostringstream& oss)
{
    auto name = comp.label();
    auto unit = comp.get_display_unit();
    oss << "start:" << name << "(" << unit << "), ";
    oss << "end:"   << name << "(" << unit << "), ";
    oss << "delta:" << name << "(" << unit << "), ";
}

// ----------------------------------------------------------------------------
// don't do anything for the runtime user bundle
template<>
void write_names(bundle_t comp, std::ostringstream& oss)
{
    (void) comp;
    (void) oss;
}

// ----------------------------------------------------------------------------
// print static component data onto the string
template<typename Tp>
void write_data(Tp comp, std::ostream& str)
{
    str << comp.get() << ", " ;
}

// ----------------------------------------------------------------------------
// don't do anything for the runtime user bundle
template<>
void write_data(bundle_t comp, std::ostream& str)
{
    (void) comp;
    (void) str;
}

// ----------------------------------------------------------------------------
// configure user global bundle
void configure_bundle()
{
    // reset any previous configuration
    bundle_t::reset();

    // get tools from environment variable
    auto env_tool = tim::get_env<std::string>("TECA_PROFILER_COMPONENTS", "");
    auto env_enum = tim::enumerate_components(tim::delimit(env_tool));
    env_enum.erase(std::remove_if(env_enum.begin(), env_enum.end(),
                                  [](int c) { return c == WALL_CLOCK; }),
                                  env_enum.end());
    // configure the bundle
    tim::configure<bundle_t>(env_enum);

    // set bundle configured
    is_configured = true;
}

} //namespace time
} //namespace measure

#endif // TECA_HAS_TIMEMORY

namespace impl
{
#if defined(TECA_ENABLE_PROFILER)

// container for data captured in a timing event
struct event
{
    event(std::string evt_name);

    // serializes the event in CSV format into the stream.
    void to_stream(std::ostream &str) const;

    enum { START=0, END=1, DELTA=2 }; // record fields

    // user provided identifier for the record
    std::string name;

    // event duration, initially start time, end time and duration
    // are recorded. when summarizing this contains min,max and sum
    // of the summariezed set of events
#if defined(TECA_HAS_TIMEMORY)
    std::vector<measurement_t> evts;
#else
    double time[3];
#endif // TECA_HAS_TIMEMORY

    // how deep is the event stack
    int depth;

    // the thread id that generated the event
    std::thread::id tid;
};

#if !defined(TECA_HAS_MPI)
using MPI_Comm = void*;
#define MPI_COMM_NULL nullptr
#endif
static MPI_Comm comm = MPI_COMM_NULL;

static int logging_enabled = 0x00;

static std::string timer_log_file = "timer.csv";

using event_log_type = std::list<impl::event>;
using thread_map_type = std::unordered_map<std::thread::id, event_log_type>;

static event_log_type event_log;
static thread_map_type active_events;
static std::mutex event_log_mutex;

// memory profiler
static teca_memory_profiler mem_prof;

#if !defined(TECA_HAS_TIMEMORY)
// return high res system time relative to system epoch
static double get_system_time()
{
    struct timeval tv;
    gettimeofday(&tv, nullptr);
    return tv.tv_sec + tv.tv_usec/1.0e6;
}
#endif // TECA_HAS_TIMEMORY

// --------------------------------------------------------------------------
event::event(std::string evt_name): name(evt_name), depth(0), tid(std::this_thread::get_id())
{
#if defined (TECA_HAS_TIMEMORY)
    if (!measure::time::is_configured)
        measure::time::configure_bundle();

    evts = { measurement_t (std::string(evt_name + "_start")), 
             measurement_t (std::string(evt_name + "_end")), 
             measurement_t (std::string(evt_name + "_delta")),
           };
#endif // TECA_HAS_TIMEMORY
}

//-----------------------------------------------------------------------------
void event::to_stream(std::ostream &str) const
{
#if defined(TECA_ENABLE_PROFILER)
    int rank = 0;
#if defined(TECA_HAS_MPI)
    int ini = 0, fin = 0;
    MPI_Initialized(&ini);
    MPI_Finalized(&fin);
    if (ini && !fin)
        MPI_Comm_rank(impl::comm, &rank);
#endif

#if !defined(TECA_HAS_TIMEMORY)
    str << rank << ", " << this->tid << ", \"" << this->name << "\", "
        << this->time[START] << ", " << this->time[END] << ", "
        << this->time[DELTA] << ", " << this->depth << std::endl;
#else
    // print all component values
    str << rank << ", " << this->tid << ", \"" << this->name << "\", ";

    // get time data and print
    auto comp_data = [&str](auto comps) { measure::time::write_data(comps, str); };

    for (auto _evnt_tuple: this->evts)
    {
        measure::apply_tuple(_evnt_tuple, comp_data);
    }

    str << this->depth << std::endl;
#endif // TECA_HAS_TIMEMORY

#else
    (void)str;
#endif // TECA_ENABLE_PROFILER
}
#endif // TECA_ENABLE_PROFILER
}




// ----------------------------------------------------------------------------
void teca_profiler::set_communicator(MPI_Comm comm)
{
#if defined(TECA_ENABLE_PROFILER) && defined(TECA_HAS_MPI)
    int ok = 0;
    MPI_Initialized(&ok);
    if (ok)
    {
        if (impl::comm != MPI_COMM_NULL)
            MPI_Comm_free(&impl::comm);

        MPI_Comm_dup(comm, &impl::comm);
    }
#else
    (void)comm;
#endif
}

// ----------------------------------------------------------------------------
void teca_profiler::set_timer_log_file(const std::string &file)
{
#if defined(TECA_ENABLE_PROFILER)
    impl::timer_log_file = file;
#else
    (void)file;
#endif
}

// ----------------------------------------------------------------------------
void teca_profiler::set_mem_prof_log_file(const std::string &file)
{
#if defined(TECA_ENABLE_PROFILER)
    impl::mem_prof.set_filename(file);
#else
    (void)file;
#endif
}

// ----------------------------------------------------------------------------
void teca_profiler::set_mem_prof_interval(int interval)
{
#if defined(TECA_ENABLE_PROFILER)
    impl::mem_prof.set_interval(interval);
#else
    (void)interval;
#endif
}

// ----------------------------------------------------------------------------
int teca_profiler::validate()
{
    int ierr = 0;
#if defined(TECA_ENABLE_PROFILER)
    if (impl::logging_enabled & 0x01)
    {
#if !defined(NDEBUG)
        impl::thread_map_type::iterator tmit = impl::active_events.begin();
        impl::thread_map_type::iterator tmend = impl::active_events.end();
        for (; tmit != tmend; ++tmit)
        {
            unsigned int n_left = tmit->second.size();
            if (n_left > 0)
            {
                std::ostringstream oss;
                impl::event_log_type::iterator it = tmit->second.begin();
                impl::event_log_type::iterator end = tmit->second.end();
                for (; it != end; ++it)
                    it->to_stream(oss);
                TECA_ERROR("Thread " << tmit->first << " has " << n_left
                    << " unmatched active events. " << std::endl
                    << oss.str())
                ierr += 1;
            }
        }
#endif
    }
#endif
    return ierr;
}

// ----------------------------------------------------------------------------
int teca_profiler::to_stream(std::ostream &os)
{
#if defined(TECA_ENABLE_PROFILER)
    if (impl::logging_enabled & 0x01)
    {
        // serialize the logged events in CSV format
        os.precision(std::numeric_limits<double>::digits10 + 2);
        os.setf(std::ios::scientific, std::ios::floatfield);

        // not locking this as it's intended to be accessed only from the main
        // thread, and all other threads are required to be finished by now
        impl::event_log_type::iterator iter = impl::event_log.begin();
        impl::event_log_type::iterator end = impl::event_log.end();

        for (; iter != end; ++iter)
            iter->to_stream(os);
    }
#else
    (void)os;
#endif
    return 0;
}

// ----------------------------------------------------------------------------
int teca_profiler::initialize()
{
#if defined(TECA_ENABLE_PROFILER)

    int rank = 0;
#if defined(TECA_HAS_MPI)
    int ok = 0;
    MPI_Initialized(&ok);
    if (ok)
    {
        // always use isolated comm space
        if (impl::comm == MPI_COMM_NULL)
            teca_profiler::set_communicator(MPI_COMM_WORLD);

        impl::mem_prof.set_communicator(impl::comm);

        MPI_Comm_rank(impl::comm, &rank);
    }
#endif

    // look for overrides in the environment
    char *tmp = nullptr;
    if ((tmp = getenv("PROFILER_ENABLE")))
        impl::logging_enabled = atoi(tmp);

    if ((tmp = getenv("PROFILER_LOG_FILE")))
        impl::timer_log_file = tmp;

    if ((tmp = getenv("MEMPROF_LOG_FILE")))
        impl::mem_prof.set_filename(tmp);

    if ((tmp = getenv("MEMPROF_INTERVAL")))
        impl::mem_prof.set_interval(atof(tmp));

    if (impl::logging_enabled & 0x02)
        impl::mem_prof.initialize();

    // report what options are in use
    if ((rank == 0) && impl::logging_enabled)
        std::cerr << "Profiler configured with event logging "
            << (impl::logging_enabled & 0x01 ? "enabled" : "disabled")
            << " and memory logging " << (impl::logging_enabled & 0x02 ? "enabled" : "disabled")
            << ", timer log file \"" << impl::timer_log_file
            << "\", memory profiler log file \"" << impl::mem_prof.get_filename()
            << "\", sampling interval " << impl::mem_prof.get_interval()
            << " seconds" << std::endl;
#endif
    return 0;
}

// ----------------------------------------------------------------------------
int teca_profiler::write_mpi_io(MPI_Comm comm, const char *file_name,
    const std::string &str)
{
#if defined(TECA_ENABLE_PROFILER) && defined(TECA_HAS_MPI)
    if (impl::logging_enabled & 0x01)
    {
        // compute the file offset
        long n_bytes = str.size();

        int rank = 0;
        int n_ranks = 1;

        MPI_Comm_rank(comm, &rank);
        MPI_Comm_size(comm, &n_ranks);

        std::vector<long> gsizes(n_ranks);
        gsizes[rank] = n_bytes;

        MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
            gsizes.data(), 1, MPI_LONG, impl::comm);

        long offset = 0;
        for (int i = 0; i < rank; ++i)
            offset += gsizes[i];

        long file_size = 0;
        for (int i = 0; i < n_ranks; ++i)
            file_size += gsizes[i];

        // write the buffer
        MPI_File fh;
        MPI_File_open(comm, file_name,
            MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

        MPI_File_set_view(fh, offset, MPI_BYTE, MPI_BYTE,
            "native", MPI_INFO_NULL);

        MPI_File_write(fh, str.c_str(), n_bytes,
            MPI_BYTE, MPI_STATUS_IGNORE);

        MPI_File_set_size(fh, file_size);

        MPI_File_close(&fh);
    }
    return 0;
#else
    (void)comm;
    (void)file_name;
    (void)str;
#endif
    return -1;
}

// ----------------------------------------------------------------------------
int teca_profiler::flush()
{
#if defined(TECA_ENABLE_PROFILER)
    std::ostringstream oss;
    teca_profiler::to_stream(oss);
    teca_profiler::write_c_stdio(impl::timer_log_file.c_str(), "a", oss.str());
    teca_profiler::validate();
    impl::event_log.clear();
#endif
    return 0;
}

// ----------------------------------------------------------------------------
int teca_profiler::write_c_stdio(const char *file_name, const char *mode,
    const std::string &str)
{
#if defined(TECA_ENABLE_PROFILER)
    if (impl::logging_enabled & 0x01)
    {
        FILE *fh = fopen(file_name, mode);
        if (!fh)
        {
            const char *estr = strerror(errno);
            TECA_ERROR("Failed to open \""
                << file_name << "\" " << estr)
            return -1;
        }

        long n_bytes = str.size();
        long nwritten = fwrite(str.c_str(), 1, n_bytes, fh);
        if (nwritten != n_bytes)
        {
            const char *estr = strerror(errno);
            TECA_ERROR("Failed to write " << n_bytes << " bytes. " << estr)
            return -1;
        }

        fclose(fh);
    }
#else
    (void)file_name;
    (void)mode;
    (void)str;
#endif
    return 0;
}

// ----------------------------------------------------------------------------
int teca_profiler::finalize()
{
#if defined(TECA_ENABLE_PROFILER)
    int ok = 0;
#if defined(TECA_HAS_MPI)
    MPI_Initialized(&ok);
#endif

    if (impl::logging_enabled & 0x01)
    {
        int rank = 0;
#if defined(TECA_HAS_MPI)
        if (ok)
            MPI_Comm_rank(impl::comm, &rank);
#endif

        // serialize the logged events in CSV format
        std::ostringstream oss;

        if (rank == 0)
#if !defined(TECA_HAS_TIMEMORY)
            oss << "# rank, thread, name, start time, end time, delta, depth" << std::endl;
#else
        {
            oss << "# rank, thread, name, ";

            measurement_t _dummy("names");
            auto comp_names = [&oss](auto comp){ measure::time::write_names(comp, oss); };
            measure::apply_tuple(_dummy, comp_names);
            oss << "depth" << std::endl;
        }
#endif // TECA_HAS_TIMEMORY

        teca_profiler::to_stream(oss);

        // free up resources
        impl::event_log.clear();

        if (ok)
            teca_profiler::write_mpi_io(impl::comm, impl::timer_log_file.c_str(), oss.str());
        else
            teca_profiler::write_c_stdio(impl::timer_log_file.c_str(), "w", oss.str());
    }

    // output the memory use profile and clean up resources
    if (impl::logging_enabled & 0x02)
        impl::mem_prof.finalize();

    // free up other resources
#if defined(TECA_HAS_MPI)
    if (ok)
        MPI_Comm_free(&impl::comm);
#endif
#endif
    return 0;
}

//-----------------------------------------------------------------------------
bool teca_profiler::enabled()
{
#if defined(TECA_ENABLE_PROFILER)
    std::lock_guard<std::mutex> lock(impl::event_log_mutex);
    return impl::logging_enabled & 0x01;
#else
    return false;
#endif
}

//-----------------------------------------------------------------------------
void teca_profiler::enable(int arg)
{
#if defined(TECA_ENABLE_PROFILER)
    std::lock_guard<std::mutex> lock(impl::event_log_mutex);
    impl::logging_enabled = arg;
#else
    (void)arg;
#endif
}

//-----------------------------------------------------------------------------
void teca_profiler::disable()
{
#if defined(TECA_ENABLE_PROFILER)
    std::lock_guard<std::mutex> lock(impl::event_log_mutex);
    impl::logging_enabled = 0x00;
#endif
}

//-----------------------------------------------------------------------------
int teca_profiler::start_event(const char* eventname)
{
#if defined(TECA_ENABLE_PROFILER)
    if (impl::logging_enabled & 0x01)
    {
        impl::event evt(eventname);
#if defined (TECA_HAS_TIMEMORY)
        evt.evts[impl::event::START].record();
        evt.evts[impl::event::DELTA].start();
#else
        evt.time[impl::event::START] = impl::get_system_time();
#endif // TECA_HAS_TIMEMORY

        std::lock_guard<std::mutex> lock(impl::event_log_mutex);
        impl::active_events[evt.tid].push_back(evt);
    }
#else
    (void)eventname;
#endif
    return 0;
}

//-----------------------------------------------------------------------------
int teca_profiler::end_event(const char* eventname)
{
#if defined(TECA_ENABLE_PROFILER)
    if (impl::logging_enabled & 0x01)
    {
#if !defined (TECA_HAS_TIMEMORY)
        // get end time
        double end_time = impl::get_system_time();
#endif // TECA_HAS_TIMEMORY

        // get this thread's event log
        std::thread::id tid = std::this_thread::get_id();

        std::lock_guard<std::mutex> lock(impl::event_log_mutex);
        impl::thread_map_type::iterator iter = impl::active_events.find(tid);
        if (iter == impl::active_events.end())
        {
            TECA_ERROR("failed to end event \"" << eventname
                << "\" thread  " << tid << " has no events")
            return -1;
        }

        impl::event evt(std::move(iter->second.back()));
        iter->second.pop_back();

#ifdef NDEBUG
        (void)eventname;
#else
        if (strcmp(eventname, evt.name.c_str()) != 0)
        {
            TECA_ERROR("Mismatched start_event/end_event. Expecting: '"
                << evt.name.c_str() << "' Got: '" << eventname << "'")
            abort();
        }
#endif
#if defined (TECA_HAS_TIMEMORY)
        evt.evts[impl::event::DELTA].stop();
        evt.evts[impl::event::END].record();
#else
        evt.time[impl::event::END] = end_time;
        evt.time[impl::event::DELTA] = end_time - evt.time[impl::event::START];
#endif // TECA_HAS_TIMEMORY
        evt.depth = iter->second.size();

        impl::event_log.emplace_back(std::move(evt));
    }
#else
    (void)eventname;
#endif
    return 0;
}
