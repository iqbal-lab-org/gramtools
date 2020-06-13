#include <boost/timer/timer.hpp>

#ifndef GRAMTOOLS_TIMER_REPORT_HPP
#define GRAMTOOLS_TIMER_REPORT_HPP

namespace gram {
class TimerReport {
 public:
  void start(std::string note);

  void stop();

  void report() const;

  template <typename TypeCol1, typename TypeCol2>
  void cout_row(TypeCol1 col1, TypeCol2 col2) const;

 private:
  using Note = std::string;
  using Entry = std::pair<Note, double>;

  Note note;
  std::vector<Entry> logger;
  boost::timer::cpu_timer timer;
};
}  // namespace gram

#endif  // GRAMTOOLS_TIMER_REPORT_HPP
