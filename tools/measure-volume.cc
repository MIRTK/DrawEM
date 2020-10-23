/*
 * Developing brain Region Annotation With Expectation-Maximization (Draw-EM)
 *
 * Copyright 2013-2016 Imperial College London
 * Copyright 2013-2016 Antonios Makropoulos
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "mirtk/Common.h"
#include "mirtk/Options.h"

#include "mirtk/IOConfig.h"
#include "mirtk/GenericImage.h"

using namespace mirtk;


// =============================================================================
// Help
// =============================================================================

// -----------------------------------------------------------------------------
void PrintHelp(const char *name)
{
	cout << endl;
	cout << "Usage: " << name << " <input>" << endl;
	cout << endl;
	cout << "Description:" << endl;
 	cout << "  Measures the volume of each label in the input image" << endl;
	cout << endl;
	cout << "Optional arguments:" << endl;
	cout << "  -voxels    Count the number of voxels instead. (default: off)" << endl;
	PrintStandardOptions(cout);
	cout << endl;
}

// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------
int main(int argc, char **argv)
{
  EXPECTS_POSARGS(1);

  InitializeIOLibrary();
  GreyImage labels(POSARG(1));

  bool voxel_count = false;
  for (ALL_OPTIONS) {
    HANDLE_BOOLEAN_OPTION("voxels", voxel_count);
    else HANDLE_STANDARD_OR_UNKNOWN_OPTION();
  }

  UnorderedMap<int, int> hist;
  const GreyPixel *label = labels.Data();
  for (int vox = 0; vox < labels.NumberOfVoxels(); ++vox, ++label) {
    if (*label != 0) {
      auto bin = hist.find(static_cast<int>(*label));
      if (bin == hist.end()) hist[*label] = 1;
      else bin->second += 1;
    }
  }

  if (voxel_count) {
    for (const auto &bin : hist) {
      cout << bin.first << " " << bin.second << "\n";
    }
  } else {
    const double vol = labels.XSize() * labels.YSize() * labels.ZSize();
    for (const auto &bin : hist) {
      cout << bin.first << " " << setprecision(8) << bin.second * vol << "\n";
    }
  }
  cout.flush();

  return 0;
}
