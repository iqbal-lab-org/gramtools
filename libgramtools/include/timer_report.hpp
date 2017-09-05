#include <boost/timer/timer.hpp>


#ifndef GRAMTOOLS_TIMER_REPORT_HPP
#define GRAMTOOLS_TIMER_REPORT_HPP

class TimerReport{
public:
    void record(std::string note);
    void report() const;

    template <typename TypeCol1, typename TypeCol2>
    void cout_row(TypeCol1 col1, TypeCol2 col2) const;

private:
    using Entry = std::pair<std::string, double>;
    std::vector<Entry> logger;
    boost::timer::cpu_timer timer;
};

#endif //GRAMTOOLS_TIMER_REPORT_HPP
