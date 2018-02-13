// The libMesh Finite Element Library.
// Copyright (C) 2002-2017 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA



// Local Includes
#include "libmesh/linear_partitioner.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/elem.h"

namespace libMesh
{

void LinearPartitioner::partition_range(MeshBase & /*mesh*/,
                                        MeshBase::element_iterator it,
                                        MeshBase::element_iterator end,
                                        const unsigned int n)
{
  libmesh_assert_greater (n, 0);

  // Check for an easy return
  if (n == 1)
    {
      this->single_partition_range (it, end);
      return;
    }

  // Create a simple linear partitioning
  LOG_SCOPE ("partition_range()", "LinearPartitioner");

  const dof_id_type blksize = std::distance(it, end) / n;

  dof_id_type e = 0;
  for ( ; it != end; ++it)
    {
      Elem * elem = *it;
      if ((e/blksize) < n)
        elem->processor_id() = cast_int<processor_id_type>(e/blksize);

      else
        elem->processor_id() = 0;

      e++;
    }
}



void LinearPartitioner::_do_partition (MeshBase & mesh,
                                       const unsigned int n)
{
  this->partition_range(mesh,
                        mesh.active_elements_begin(),
                        mesh.active_elements_end(),
                        n);
}

} // namespace libMesh
