/**
 * \file ColorCoder.h \brief
 *
 * Copyright 2007-2022 IMP Inventors. All rights reserved.
 *
 */

#ifndef FOXS_COLOR_CODER_H
#define FOXS_COLOR_CODER_H

#include "foxs_config.h"
#include <cstdio>
#include <stdio.h>

namespace foxs { namespace internal {

class ColorCoder {
 public:
  static void get_color_for_id(int &r, int &g, int &b, int id);

  static void set_number(int number) {
    if (number < 6) {
      diff_ = 150;
      return;
    }
    diff_ = 150 / (number / 6);
  }

  static void html_hex_color(char *out_color, int id) {
    int r, g, b;
    get_color_for_id(r, g, b, id);
    sprintf(out_color, "%02X%02X%02X", r, g, b);
  }

  static void jmol_dec_color(char *out_color, int id) {
    int r, g, b;
    get_color_for_id(r, g, b, id);
    sprintf(out_color, "[%d ,%d ,%d]", r, g, b);
  }

  static int diff_;
};

} }  // namespace foxs::internal

#endif /* FOXS_COLOR_CODER_H */
