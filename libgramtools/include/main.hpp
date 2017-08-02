#ifndef GRAMTOOLS_MAIN_HPP
#define GRAMTOOLS_MAIN_HPP


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

Parameters parse_command_line_parameters(int argc, const char *const *argv);


#endif //GRAMTOOLS_MAIN_HPP
