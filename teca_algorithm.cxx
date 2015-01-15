#include "teca_algorithm.h"
#include "teca_dataset.h"
#include "teca_algorithm_executive.h"
#include "teca_threadsafe_queue.h"

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <utility>
#include <algorithm>
#include <thread>
#include <mutex>
#include <atomic>

using std::vector;
using std::map;
using std::string;
using std::ostream;
using std::istream;
using std::pair;
using std::thread;
using std::mutex;
using std::lock_guard;
using std::atomic;

namespace {
// convenience functions for accessing port and algorithm
// from an output port
static
p_teca_algorithm &algorithm(teca_algorithm::output_port_t &op)
{ return op.first; }

static
unsigned int &port(teca_algorithm::output_port_t &op)
{ return op.second; }
};


class teca_algorithm_thread_pool;
typedef std::shared_ptr<teca_algorithm_thread_pool> p_teca_algorithm_thread_pool;

// a class to manage a fixed size pool of threads that dispatch
// data requests to teca_algorithm
class teca_algorithm_thread_pool
{
public:
    typedef
    pair<teca_algorithm::output_port_t, teca_meta_data>
        data_request_t;

    // construct the thread pool with the default number
    // of threads (1-the number of cores)
    teca_algorithm_thread_pool();
    // contruct the thread pool with n worker threads
    teca_algorithm_thread_pool(unsigned int n);

    ~teca_algorithm_thread_pool();

    // add a task
    void add_data_request(
        teca_algorithm::output_port_t &port,
        const teca_meta_data &req);

    // wait for all pending tasks to complete
    void wait_pending();

    // get the number of threads
    unsigned int size()
    { return m_threads.size(); }

protected:
    teca_algorithm_thread_pool(const teca_algorithm_thread_pool &) = delete;
    void operator=(const teca_algorithm_thread_pool &) = delete;

    // create n threads for the pool
    void create_threads(unsigned int n_threads);

private:
    atomic<bool> m_live;
    atomic<unsigned long> m_pending;
    teca_threadsafe_queue<data_request_t> m_queue;
    vector<thread> m_threads;
};

namespace {
// convenience functions for accessing output port and request
// from a task meta data
teca_algorithm::output_port_t &port(teca_algorithm_thread_pool::data_request_t &task)
{ return task.first; }

teca_meta_data &request(teca_algorithm_thread_pool::data_request_t &task)
{ return task.second; }
};

// --------------------------------------------------------------------------
teca_algorithm_thread_pool::teca_algorithm_thread_pool()
    : m_live(true), m_pending(0)
{
    unsigned int n_threads = std::max(1u, thread::hardware_concurrency()-1);
    this->create_threads(n_threads);
}

// --------------------------------------------------------------------------
teca_algorithm_thread_pool::teca_algorithm_thread_pool(unsigned int n_threads)
    : m_live(true), m_pending(0)
{
    this->create_threads(n_threads);
}

// --------------------------------------------------------------------------
void teca_algorithm_thread_pool::create_threads(unsigned int n_threads)
{
    for (unsigned int i = 0; i < n_threads; ++i)
    {
        m_threads.push_back(thread([this]()
        {
            while (m_live.load())
            {
                data_request_t task;
                if (m_queue.try_pop(task))
                {
                    teca_algorithm::request_data(
                        ::port(task), ::request(task));
                    --m_pending;
                }
                else
                {
                    std::this_thread::yield();
                }
            }
        }));
    }
}

// --------------------------------------------------------------------------
teca_algorithm_thread_pool::~teca_algorithm_thread_pool()
{
    m_live = false;
    std::for_each(
        m_threads.begin(), m_threads.end(),
        [](thread &t) { t.join(); });
}

// --------------------------------------------------------------------------
void teca_algorithm_thread_pool::add_data_request(
        teca_algorithm::output_port_t &port,
        const teca_meta_data &req)
{
    m_queue.push(data_request_t(port, req));
    ++m_pending;
}

// --------------------------------------------------------------------------
void teca_algorithm_thread_pool::wait_pending()
{
    while (m_pending.load())
        std::this_thread::yield();
}






// implementation for managing input connections
// and cached output data
class teca_algorithm_internals
{
public:
    teca_algorithm_internals();
    ~teca_algorithm_internals();

    teca_algorithm_internals(const teca_algorithm_internals &other) = delete;
    teca_algorithm_internals(teca_algorithm_internals &&other) = delete;

    teca_algorithm_internals &operator=(const teca_algorithm_internals &other) = delete;
    teca_algorithm_internals &operator=(teca_algorithm_internals &&other) = delete;

    // this setsup the output cache. calling
    // this will remove all cached data.
    void set_number_of_inputs(unsigned int n);
    void set_number_of_outputs(unsigned int n);

    unsigned int get_number_of_inputs() const;
    unsigned int get_number_of_outputs() const;

    // set/get the input
    teca_algorithm::output_port_t &get_input(unsigned int i);

    void set_input(
        unsigned int conn,
        const teca_algorithm::output_port_t &port);

    // insert a dataset into the cache for the
    // given request (thread safe)
    int cache_output_data(
        unsigned int port,
        const teca_meta_data &request,
        p_teca_dataset data);

    // get a pointer to the cached dataset. if
    // a no dataset is cached for the given request
    // then return is null (thread safe)
    p_teca_dataset get_output_data(
        unsigned int port,
        const teca_meta_data &request);

    // get a pointer to the "first" cached dataset.
    p_teca_dataset get_output_data(unsigned int port);

    // clear the cache (thread safe)
    void clear_data_cache(unsigned int port);

    // remove a dataset at the top or bottom of the cache
    void pop_cache(unsigned int port, int top);

    // sets the maximum nuber of datasets to cache
    // per output port
    void set_data_cache_size(unsigned int n);
    unsigned int get_data_cache_size() const;

    // set/clear modified flag for the given port
    void set_modified();
    void set_modified(unsigned int port);

    void clear_modified();
    void clear_modified(unsigned int port);

    // get the modified state of the given port
    int get_modified(unsigned int port) const;

    // set/get the executive
    p_teca_algorithm_executive get_executive();
    void set_executive(p_teca_algorithm_executive &exec);

    // set the number of threads
    void set_number_of_threads(unsigned int n_threads);

    // print internal state to the stream
    void to_stream(ostream &os) const;
    void from_stream(istream &is);

    // algorithm description
    string name;

    // links to upstream stages indexed by input
    // connection id
    vector<teca_algorithm::output_port_t> inputs;

    // cached output data. maps from a request
    // to the cached dataset, one per output port
    unsigned int data_cache_size;
    typedef map<teca_meta_data, p_teca_dataset> req_data_map;
    vector<req_data_map> data_cache;
    mutex data_cache_mutex;

    // flag that indicates if the cache on output port
    // i is invalid
    vector<int> modified;

    // executive
    p_teca_algorithm_executive exec;

    // thread pool
    p_teca_algorithm_thread_pool thread_pool;
};

// --------------------------------------------------------------------------
teca_algorithm_internals::teca_algorithm_internals()
            :
    name("teca_algorithm"),
    data_cache_size(1),
    modified(1),
    exec(teca_algorithm_executive::New()),
    thread_pool(new teca_algorithm_thread_pool(1))
{
    this->set_number_of_outputs(1);
}

// --------------------------------------------------------------------------
teca_algorithm_internals::~teca_algorithm_internals()
{}
/*
// --------------------------------------------------------------------------
teca_algorithm_internals::teca_algorithm_internals(
    const teca_algorithm_internals &other)
            :
    name(other.name),
    inputs(other.inputs),
    data_cache_size(other.data_cache_size),
    data_cache(other.data_cache),
    modified(other.modified),
    exec(other.exec),
    request_queue(other.request_queue)
{}

// --------------------------------------------------------------------------
teca_algorithm_internals::teca_algorithm_internals(
    teca_algorithm_internals &&other)
{
    this->name.swap(other.name);
    this->inputs.swap(other.inputs);
    this->data_cache_size = other.data_cache_size;
    this->data_cache.swap(other.data_cache);
    this->modified.swap(other.modified);
    this->exec.swap(other.exec);
    this->equest_queue.swap(other.request_queue);
}

// --------------------------------------------------------------------------
teca_algorithm_internals &teca_algorithm_internals::operator=(
    const teca_algorithm_internals &other)
{
    if (&other == this)
        return *this;

    this->name = other.name;
    this->inputs = other.inputs;
    this->data_cache_size = other.data_cache_size;
    this->data_cache = other.data_cache;
    this->modified = other.modified;
    this->exec = other.exec;
    this->request_queue = other.request_queue;

    return *this;
}

// --------------------------------------------------------------------------
teca_algorithm_internals &teca_algorithm_internals::operator=(
    teca_algorithm_internals &&other)
{
    if (&other == this)
        return *this;

    this->name.clear();
    this->name.swap(other.name);

    this->inputs.clear();
    this->inputs.swap(other.inputs);

    this->data_cache_size = other.data_cache_size;

    this->data_cache.clear();
    this->data_cache.swap(other.data_cache);

    this->modified.clear();
    this->modified.swap(other.modified);

    this->exec.reset();
    this->exec.swap(other.exec);

    this->request_queue.clear();
    this->request_queue.swap(other.request_queue);

    return *this;
}
*/

// --------------------------------------------------------------------------
void teca_algorithm_internals::set_number_of_inputs(unsigned int n)
{
    this->inputs.clear();
    this->inputs.resize(n, teca_algorithm::output_port_t(nullptr, 0));
}

// --------------------------------------------------------------------------
unsigned int teca_algorithm_internals::get_number_of_inputs() const
{
    return this->inputs.size();
}

// --------------------------------------------------------------------------
void teca_algorithm_internals::set_input(
    unsigned int conn,
    const teca_algorithm::output_port_t &port)
{
    this->inputs[conn] = port;
}

// --------------------------------------------------------------------------
teca_algorithm::output_port_t &teca_algorithm_internals::get_input(
    unsigned int conn)
{
    return this->inputs[conn];
}

// --------------------------------------------------------------------------
void teca_algorithm_internals::set_number_of_outputs(unsigned int n)
{
    if (n < 1)
    {
        TECA_ERROR("invalid number of outputs " << n)
        n = 1;
    }

    // create a chacne for each output
    this->data_cache.clear();
    this->data_cache.resize(n);

    // create a modified flag for each output
    this->modified.clear();
    this->modified.resize(n, 1);
}

// --------------------------------------------------------------------------
unsigned int teca_algorithm_internals::get_number_of_outputs() const
{
    return this->data_cache.size();
}

// --------------------------------------------------------------------------
void teca_algorithm_internals::set_data_cache_size(unsigned int n)
{
    this->data_cache_size = n;
    unsigned int n_out = this->get_number_of_outputs();
    for (unsigned int i = 0; i < n_out; ++i)
    {
        req_data_map &cache = this->data_cache[i];

        while (cache.size() >  n)
            cache.erase(cache.begin());
    }
}

// --------------------------------------------------------------------------
unsigned int teca_algorithm_internals::get_data_cache_size() const
{
    return this->data_cache_size;
}

// --------------------------------------------------------------------------
void teca_algorithm_internals::clear_data_cache(unsigned int port)
{
    std::lock_guard<mutex> lock(this->data_cache_mutex);

    this->data_cache[port].clear();
}

// --------------------------------------------------------------------------
void teca_algorithm_internals::pop_cache(unsigned int port, int top)
{
    req_data_map &cache = this->data_cache[port];

    if (cache.empty())
        return;

    if (top)
        cache.erase(--cache.end()); // newest
    else
        cache.erase(cache.begin()); // oldest
}

// --------------------------------------------------------------------------
int teca_algorithm_internals::cache_output_data(
    unsigned int port,
    const teca_meta_data &request,
    p_teca_dataset data)
{
    if (this->data_cache_size)
    {
        std::lock_guard<mutex> lock(this->data_cache_mutex);

        req_data_map &cache = this->data_cache[port];

        auto res = cache.insert(req_data_map::value_type(request, data));
        if (!res.second)
            res.first->second = data;

        while (cache.size() >= this->data_cache_size)
            cache.erase(cache.begin());
    }
    return 0;
}

// --------------------------------------------------------------------------
p_teca_dataset teca_algorithm_internals::get_output_data(
    unsigned int port,
    const teca_meta_data &request)
{
    std::lock_guard<mutex> lock(this->data_cache_mutex);

    req_data_map &cache = this->data_cache[port];

    req_data_map::iterator it = cache.find(request);
    if (it != cache.end())
    {
        return it->second;
    }

    return p_teca_dataset();
}

// --------------------------------------------------------------------------
p_teca_dataset teca_algorithm_internals::get_output_data(unsigned int port)
{
    std::lock_guard<mutex> lock(this->data_cache_mutex);

    req_data_map &cache = this->data_cache[port];

    if (cache.empty())
        return p_teca_dataset();

    return cache.rbegin()->second; // newest
}

// --------------------------------------------------------------------------
int teca_algorithm_internals::get_modified(unsigned int port) const
{
    return this->modified[port];
}

// --------------------------------------------------------------------------
void teca_algorithm_internals::set_modified(unsigned int port)
{
    this->modified[port] = 1;
}

// --------------------------------------------------------------------------
void teca_algorithm_internals::set_modified()
{
    std::for_each(
        this->modified.begin(), this->modified.end(),
        [](int &m){ m = 1; });
}

// --------------------------------------------------------------------------
void teca_algorithm_internals::clear_modified(unsigned int port)
{
    this->modified[port] = 0;
}

// --------------------------------------------------------------------------
void teca_algorithm_internals::clear_modified()
{
    std::for_each(
        this->modified.begin(), this->modified.end(),
        [](int &m){ m = 0; });
}

// --------------------------------------------------------------------------
p_teca_algorithm_executive teca_algorithm_internals::get_executive()
{
    return this->exec;
}

// --------------------------------------------------------------------------
void teca_algorithm_internals::set_executive(p_teca_algorithm_executive &e)
{
    this->exec = e;
}

// --------------------------------------------------------------------------
void teca_algorithm_internals::set_number_of_threads(unsigned int n_threads)
{
    if (this->thread_pool->size() != n_threads)
    {
        this->thread_pool
            = std::make_shared<teca_algorithm_thread_pool>(n_threads);
    }
}

// --------------------------------------------------------------------------
void teca_algorithm_internals::to_stream(ostream &os) const
{
    // TODO
    (void) os;
}

// --------------------------------------------------------------------------
void teca_algorithm_internals::from_stream(istream &is)
{
    // TODO
    (void) is;
}










// --------------------------------------------------------------------------
p_teca_algorithm teca_algorithm::New()
{
    return p_teca_algorithm(new teca_algorithm);
}

// --------------------------------------------------------------------------
teca_algorithm::teca_algorithm()
{
    this->internals = new teca_algorithm_internals;
}

// --------------------------------------------------------------------------
teca_algorithm::~teca_algorithm()
{
    delete this->internals;
}

/*// --------------------------------------------------------------------------
teca_algorithm::teca_algorithm(const teca_algorithm &src)
    : internals(new teca_algorithm_internals(*src.internals))
{}

// --------------------------------------------------------------------------
teca_algorithm::teca_algorithm(teca_algorithm &&src)
    : internals(new teca_algorithm_internals(std::move(*src.internals)))
{}

// --------------------------------------------------------------------------
teca_algorithm &teca_algorithm::operator=(const teca_algorithm &src)
{
    if (this == &src)
        return *this;

    teca_algorithm_internals *tmp
        = new teca_algorithm_internals(*src.internals);

    delete this->internals;
    this->internals = tmp;

    return *this;
}

// --------------------------------------------------------------------------
teca_algorithm &teca_algorithm::operator=(teca_algorithm &&src)
{
    if (this == &src)
        return *this;

    *this->internals = std::move(*src.internals);

    return *this;
}*/

// --------------------------------------------------------------------------
teca_algorithm::output_port_t teca_algorithm::get_output_port(
    unsigned int port)
{
    return teca_algorithm::output_port_t(this->shared_from_this(), port);
}

// --------------------------------------------------------------------------
void teca_algorithm::set_input_connection(
    unsigned int conn,
    const output_port_t &upstream)
{
    this->internals->set_input(conn, upstream);
    this->set_modified();
}

// --------------------------------------------------------------------------
void teca_algorithm::remove_input_connection(unsigned int id)
{
    this->set_input_connection(id, output_port_t(nullptr, 0));
    this->set_modified();
}

// --------------------------------------------------------------------------
void teca_algorithm::clear_input_connections()
{
    unsigned int n = this->internals->get_number_of_inputs();
    for (unsigned int i=0; i<n; ++i)
    {
        this->remove_input_connection(i);
    }
    this->set_modified();
}

// --------------------------------------------------------------------------
p_teca_dataset teca_algorithm::get_output_data(
    unsigned int port,
    const teca_meta_data &request,
    int filter)
{
    teca_meta_data key;
    if (filter)
        key = this->get_cache_key(port, request);
    else
        key = request;

    return this->internals->get_output_data(port, key);
}

// --------------------------------------------------------------------------
p_teca_dataset teca_algorithm::get_output_data(unsigned int port)
{
    return this->internals->get_output_data(port);
}

// --------------------------------------------------------------------------
void teca_algorithm::pop_cache(unsigned int port, int top)
{
    this->internals->pop_cache(port, top);
    this->set_modified(port);
}

// --------------------------------------------------------------------------
void teca_algorithm::set_cache_size(unsigned int n)
{
    if (n == this->internals->get_data_cache_size())
        return;

    this->internals->set_data_cache_size(n);
    this->set_modified();
}

// --------------------------------------------------------------------------
void teca_algorithm::set_executive(p_teca_algorithm_executive exec)
{
    this->internals->set_executive(exec);
}

// --------------------------------------------------------------------------
void teca_algorithm::set_number_of_inputs(unsigned int n)
{
    this->internals->set_number_of_inputs(n);
}

// --------------------------------------------------------------------------
void teca_algorithm::set_number_of_outputs(unsigned int n)
{
    this->internals->set_number_of_outputs(n);
}

// --------------------------------------------------------------------------
void teca_algorithm::set_number_of_threads(unsigned int n)
{
    this->internals->set_number_of_threads(n);
}

// --------------------------------------------------------------------------
teca_meta_data teca_algorithm::get_output_meta_data(
    unsigned int port,
    const vector<teca_meta_data> &input_md)
{
    // the default implementation passes meta data through
    if (input_md.size())
        return input_md[0];

    return teca_meta_data();
}

// --------------------------------------------------------------------------
p_teca_dataset teca_algorithm::execute(
        unsigned int port,
        const vector<p_teca_dataset> &input_data,
        const teca_meta_data &request)
{
    (void) port;
    (void) input_data;
    (void) request;

    // default implementation does nothing
    return p_teca_dataset();
}

// --------------------------------------------------------------------------
vector<teca_meta_data> teca_algorithm::get_upstream_request(
    unsigned int port,
    const vector<teca_meta_data> &input_md,
    const teca_meta_data &request)
{
    (void) port;
    (void) input_md;

    // default implementation forwards request upstream
    return vector<teca_meta_data>(
        this->internals->get_number_of_inputs(), request);
}

// --------------------------------------------------------------------------
teca_meta_data teca_algorithm::get_cache_key(
    unsigned int port,
    const teca_meta_data &request) const
{
    // default implementation passes the request through
    return request;
}

// --------------------------------------------------------------------------
teca_meta_data teca_algorithm::get_output_meta_data(output_port_t &current)
{
    p_teca_algorithm alg = ::algorithm(current);
    unsigned int port = ::port(current);

    // gather upstream metadata one per input
    // connection.

    unsigned int n_inputs = alg->internals->get_number_of_inputs();
    vector<teca_meta_data> input_md(n_inputs);
    for (unsigned int i = 0; i < n_inputs; ++i)
    {
        input_md[i]
          = alg->get_output_meta_data(alg->internals->get_input(i));
    }

    // now that we have metadata for the algorithm's
    // inputs, call the override to do the actual work
    // of reporting output meta data
    return alg->get_output_meta_data(port, input_md);
}

// --------------------------------------------------------------------------
p_teca_dataset teca_algorithm::request_data(
    output_port_t &current,
    const teca_meta_data &request)
{
    // execute current algorithm to fulfill the request.
    // return the data
    p_teca_algorithm &alg = ::algorithm(current);
    unsigned int port = ::port(current);

    // check for cached data
    teca_meta_data key = alg->get_cache_key(port, request);
    p_teca_dataset out_data = alg->internals->get_output_data(port, key);
    if (!out_data)
    {
        // determine what data is available on our inputs
        unsigned int n_inputs = alg->internals->get_number_of_inputs();
        vector<teca_meta_data> input_md(n_inputs);
        for (unsigned int i = 0; i < n_inputs; ++i)
        {
            input_md[i]
              = alg->get_output_meta_data(alg->internals->get_input(i));
        }

        // get requests for upstream data
        vector<teca_meta_data> up_reqs
            = alg->get_upstream_request(port, input_md, request);

        // get the upstream data mapping the requests round-robbin
        // on to the inputs
        size_t n_up_reqs = up_reqs.size();
        vector<p_teca_dataset> input_data(n_up_reqs);
        for (unsigned int i = 0; i < n_up_reqs; ++i)
        {
            if (!up_reqs[i].empty())
            {
                input_data[i]
                  = alg->request_data(alg->internals->get_input(i%n_inputs), up_reqs[i]);
            }
        }

        // execute override
        out_data = alg->execute(port, input_data, request);

        // cache
        alg->internals->cache_output_data(port, key, out_data);
    }

    return out_data;
}

// --------------------------------------------------------------------------
int teca_algorithm::validate_cache(output_port_t &current)
{
    p_teca_algorithm alg = ::algorithm(current);
    unsigned int port = ::port(current);

    unsigned int n = alg->internals->get_number_of_inputs();
    for (unsigned int i = 0; i<n; ++i)
    {
        output_port_t upstream = alg->internals->get_input(i);

        if (alg->validate_cache(upstream)
            || alg->internals->get_modified(port))
        {
            alg->internals->clear_data_cache(port);
            return 1;
        }
    }

    // only if no upstream has been modified can we
    // report the cache is valid
    return 0;
}

// --------------------------------------------------------------------------
void teca_algorithm::clear_modified(output_port_t current)
{
    p_teca_algorithm alg = ::algorithm(current);
    unsigned int port = ::port(current);

    unsigned int n = alg->internals->get_number_of_inputs();
    for (unsigned int i = 0; i < n; ++i)
    {
        alg->clear_modified(alg->internals->get_input(i));
    }

    alg->internals->clear_modified(::port(current));
}

// --------------------------------------------------------------------------
int teca_algorithm::update()
{
    p_teca_algorithm_executive exec = this->internals->get_executive();

    // produce data on each of our outputs
    unsigned int n_out = this->internals->get_number_of_outputs();
    for (unsigned int i = 0; i < n_out; ++i)
    {
        output_port_t port = this->get_output_port(i);

        // make sure caches are wiped where inputs have changed
        teca_algorithm::validate_cache(port);

        // initialize the executive
        if (exec->initialize(teca_algorithm::get_output_meta_data(port)))
        {
            TECA_ERROR("failed to initialize the executive")
            return -1;
        }

        // let thread pool execute the requests
        teca_meta_data request;
        while ((request = exec->get_next_request()))
            this->internals->thread_pool->add_data_request(port, request);

        // wait for execution of all requests to complete
        this->internals->thread_pool->wait_pending();

        // clear modfied flags
        teca_algorithm::clear_modified(port);
    }
    return 0;
}

// --------------------------------------------------------------------------
void teca_algorithm::set_modified()
{
    this->internals->set_modified();
}

// --------------------------------------------------------------------------
void teca_algorithm::set_modified(unsigned int port)
{
    this->internals->set_modified(port);
}

// --------------------------------------------------------------------------
void teca_algorithm::to_stream(ostream &os) const
{
    this->internals->to_stream(os);
}

// --------------------------------------------------------------------------
void teca_algorithm::from_stream(istream &is)
{
    this->internals->from_stream(is);
}
