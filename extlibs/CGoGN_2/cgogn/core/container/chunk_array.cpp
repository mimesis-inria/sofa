/*******************************************************************************
* CGoGN: Combinatorial and Geometric modeling with Generic N-dimensional Maps  *
* Copyright (C) 2015, IGG Group, ICube, University of Strasbourg, France       *
*                                                                              *
* This library is free software; you can redistribute it and/or modify it      *
* under the terms of the GNU Lesser General Public License as published by the *
* Free Software Foundation; either version 2.1 of the License, or (at your     *
* option) any later version.                                                   *
*                                                                              *
* This library is distributed in the hope that it will be useful, but WITHOUT  *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or        *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License  *
* for more details.                                                            *
*                                                                              *
* You should have received a copy of the GNU Lesser General Public License     *
* along with this library; if not, write to the Free Software Foundation,      *
* Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.           *
*                                                                              *
* Web site: http://cgogn.unistra.fr/                                           *
* Contact information: cgogn@unistra.fr                                        *
*                                                                              *
*******************************************************************************/

#define CGOGN_CORE_CONTAINER_CHUNK_ARRAY_CPP_

#include <cgogn/core/dll.h>
#include <cgogn/core/container/chunk_array.h>

namespace cgogn
{

//template class CGOGN_CORE_API ChunkArray<CGOGN_CHUNK_SIZE, bool>;
template class CGOGN_CORE_API ChunkArray<CGOGN_CHUNK_SIZE, uint32>;
template class CGOGN_CORE_API ChunkArray<CGOGN_CHUNK_SIZE, uint8>;
template class CGOGN_CORE_API ChunkArray<CGOGN_CHUNK_SIZE, std::array<float32, 3>>;
template class CGOGN_CORE_API ChunkArray<CGOGN_CHUNK_SIZE, std::array<float64, 3>>;
template class CGOGN_CORE_API ChunkArrayBool<CGOGN_CHUNK_SIZE>;

} // namespace cgogn