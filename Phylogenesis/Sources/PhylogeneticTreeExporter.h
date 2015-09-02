/*  This file is part of Victor.

    Victor is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Victor is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Victor.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef __PhylogeneticTreeExporter_H__
#define __PhylogeneticTreeExporter_H__

#include <Cluster.h>

namespace Victor { namespace Phylogenesis{

	/** @brief   Abstract class for clustering tree exporter.
	*   @author  Andrea Nasari (andrea.nasari\@studenti.unipd.it)
	*
	**/
	class PhylogeneticTreeExporter
	{
	public:
		/// Destructor.
		virtual ~PhylogeneticTreeExporter(){}

		/// Export the clustering tree
		virtual void save(Cluster* tree) = 0;
	};

}}

#endif
