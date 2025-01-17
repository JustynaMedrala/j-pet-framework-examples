/**
 *  @copyright Copyright 2018 The J-PET Framework Authors. All rights reserved.
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may find a copy of the License in the LICENCE file.
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 *  @file main.cpp
 */

#include <JPetManager/JPetManager.h>
#include "../LargeBarrelAnalysis/EventFinder.h"
#include "EventAnalyzer.h"
using namespace std;

int main(int argc, const char* argv[])
{
  try {
    JPetManager& manager = JPetManager::getManager();
    
    manager.registerTask<EventFinder>("EventFinder");
    //manager.registerTask<EventAnalyzer>("EventAnalyzer");
    
    manager.useTask("EventFinder", "hits", "unk.evt");
    //manager.useTask("EventAnalyzer", "unk.evt", "ana.evt");
    
    manager.run(argc, argv);
  } catch (const std::exception& except) {
    std::cerr << "Unrecoverable error occured:" << except.what() << "Exiting the program!" << std::endl;
    return EXIT_FAILURE;
  }
}
