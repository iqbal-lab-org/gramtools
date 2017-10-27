#include <iostream>
#include <iomanip>
#include <vector>
#include <boost/timer/timer.hpp>

#include "timer_report.hpp"


void TimerReport::start(std::string note) {
    this->note = note;
    timer.start();
}


void TimerReport::stop() {
    if (this->note.empty())
        std::cerr << "TimerReport stop called with empty note" << std::endl;
    boost::timer::cpu_times times = timer.elapsed();
    double elapsed_time = (times.user + times.system) * 1e-9;
    Entry entry = std::make_pair(note, elapsed_time);
    logger.push_back(entry);
    this->note = "";
}


void TimerReport::report() const {
    std::cout << "\nTimer report:" << std::endl;
    cout_row(" ", "seconds");

    for (const auto &entry: TimerReport::logger) {
        auto &note = std::get<0>(entry);
        auto &elapsed_time = std::get<1>(entry);
        cout_row(note, elapsed_time);
    }
}


template<typename TypeCol1, typename TypeCol2>
void TimerReport::cout_row(TypeCol1 col1, TypeCol2 col2) const {
    std::cout << std::setw(20) << std::right << col1
              << std::setw(10) << std::right << col2
              << std::endl;
}