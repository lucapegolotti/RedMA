// Reduced Modeling of Arteries (ReMA)
// Copyright (C) 2019  Luca Pegolotti
//
// ReMA is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// ReMA is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include <stdio.h>
#include <tinyxml2.h>

// Simple example that shows how to use tinyxlm2 to read XLM files

int main(int argc, char **argv)
{
    tinyxml2::XMLDocument doc;
    doc.LoadFile("data/test.xml");

    const char* title = doc.FirstChildElement("rootnode")->
                            FirstChildElement("buildingblock")->
                            FirstChildElement("x0")->GetText();
    printf( "First attribute: %s\n", title );

    const char* name = doc.FirstChildElement("rootnode")->
                            FirstChildElement()->Value();
    printf( "Name of first child of rootnode: %s\n", name );

    return 0;
}
