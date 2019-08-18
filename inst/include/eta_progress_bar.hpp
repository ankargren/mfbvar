/*
 * eta_progress_bar.hpp
 *
 * A custom ProgressBar class to display a progress bar with time estimation
 *
 * Author: clemens@nevrome.de
 *
 */
#ifndef _RcppProgress_ETA_PROGRESS_BAR_HPP
#define _RcppProgress_ETA_PROGRESS_BAR_HPP

#include <R_ext/Print.h>
#include <ctime>
#include <stdio.h>
#include <sstream>
#include <string.h>

#include "progress_bar.hpp"

// for unices only
#if !defined(WIN32) && !defined(__WIN32) && !defined(__WIN32__)
#include <Rinterface.h>
#endif

class ETAProgressBar: public ProgressBar{
  public: // ====== LIFECYCLE =====

    /**
    * Main constructor
    */
    ETAProgressBar()  {
      _max_ticks = 50;
      _finalized = false;
      _timer_flag = true;
    }

    ~ETAProgressBar() {
    }

    public: // ===== main methods =====

      void display() {
      }

      // update display
      void update(float progress) {

        // stop if already finalized
        if (_finalized) return;

          // create progress bar string
          std::string progress_bar_string = _current_ticks_display(progress);

          // ensure overwriting of old time info
          int empty_length = time_string.length();
          std::string empty_space = std::string(empty_length, ' ');

          // merge progress bar and time string
          std::stringstream strs;
          strs << "0% [" << progress_bar_string << "] " << pct <<"%";
          std::string temp_str = strs.str();
          char const* char_type = temp_str.c_str();

          // print: remove old and replace with new
          REprintf("\r");
          REprintf("%s", char_type);

          // finalize display when ready
          if(progress == 1) {
            _finalize_display();
          }
      }

      void end_display() {
        update(1);
      }

      protected: // ==== other instance methods =====

        // update the ticks display corresponding to progress
        std::string _current_ticks_display(float progress) {

          int nb_ticks = _compute_nb_ticks(progress);

          std::string cur_display = _construct_ticks_display_string(nb_ticks);

          return cur_display;
        }

        // construct progress bar display
        std::string _construct_ticks_display_string(int nb) {

          std::stringstream ticks_strs;
          for (int i = 0; i < (_max_ticks - 1); ++i) {
            if (i < nb) {
              ticks_strs << "+";
            } else {
              ticks_strs << "-";
            }
          }
          std::string tick_space_string = ticks_strs.str();

          return tick_space_string;
        }

        // finalize
        void _finalize_display() {
          if (_finalized) return;

          REprintf("\n");
          flush_console();
          _finalized = true;
        }

        // compute number of ticks according to progress
        int _compute_nb_ticks(float progress) {
          return int(progress * _max_ticks);
        }

        // N.B: does nothing on windows
        void flush_console() {
#if !defined(WIN32) && !defined(__WIN32) && !defined(__WIN32__)
          R_FlushConsole();
#endif
        }

        private: // ===== INSTANCE VARIABLES ====
          int _max_ticks;   		// the total number of ticks to print
          bool _finalized;
          bool _timer_flag;
          time_t start,end;

};

#endif
