#include "sdsl/suffix_arrays.hpp"

//need to skip chromosome separators as well

using namespace sdsl;

bool skip(const FM_Index &fm_index,
          uint64_t &left, uint64_t &right,
          uint64_t &left_rev, uint64_t &right_rev,
          uint64_t num, uint64_t maxx) {
    
    uint64_t site_end, site_start, ind_start, ind_end;
    bool last = false;

    assert(left < right);
    assert(right <= fm_index.size());
    assert(num > 4);

    if (num % 2 == 1) {
        uint64_t num_begin = fm_index.C[fm_index.char2comp[num]];

        site_end = std::max(fm_index[num_begin], fm_index[num_begin + 1]);
        site_start = std::min(fm_index[num_begin], fm_index[num_begin + 1]);

        if (fm_index[num_begin] > fm_index[num_begin + 1]) {
            ind_end = num_begin;
            ind_start = num_begin + 1;
        } else {
            ind_end = num_begin + 1;
            ind_start = num_begin;
        }

        if (right - left == 1) {
            if (fm_index[left] == site_end + 1) {
                last = true;

                left = num_begin;
                if (num + 2 <= maxx) right = fm_index.C[fm_index.char2comp[num + 2]];
                else right = fm_index.size();

                left_rev = left;
                right_rev = right;
            } else {
                //site=num;
                //allele.push_back(1);

                left = ind_start;
                right = ind_start + 1;

                left_rev = right;
                right_rev = left;
                //	first=T;
            }
        } else {
            //site=num;
            //allele.push_back(1);

            left = num_begin;
            if (num + 2 <= maxx) right = fm_index.C[fm_index.char2comp[num + 2]];
            else right = fm_index.size();

            left_rev = left;
            right_rev = right;
            //first=T;
        }
    } else {
        uint64_t num_begin = fm_index.C[fm_index.char2comp[num - 1]];

        site_start = std::min(fm_index[num_begin], fm_index[num_begin + 1]);
        site_end = std::max(fm_index[num_begin], fm_index[num_begin + 1]);

        /*
        if (right-left==1) {
          site=num-1;
          allele.push_back(mask_a[fm_index[left]]);
        }
        else {
          site=num-1;
          for (i=left;i<right;i++)
        allele.push_back(mask_a[fm_index[i]]);
        }
        */
        if (fm_index[num_begin] < fm_index[num_begin + 1]) {
            left = num_begin;
            right = num_begin + 1;

            //left_rev=site_end;
            //right_rev=site_end+1;
            left_rev = num_begin + 1;
            right_rev = num_begin + 2;
        } else {
            left = num_begin + 1;
            right = num_begin + 2;

            //left_rev=site_end;
            //right_rev=site_end+1;
            left_rev = num_begin;
            right_rev = num_begin + 1;
        }
    }

    assert(right > left || !(std::cerr << "False: " << right << ">" << left << " " << num << " " << last << " "));
    //  assert(right_rev-left_rev == right-left);

    return (last);
}
