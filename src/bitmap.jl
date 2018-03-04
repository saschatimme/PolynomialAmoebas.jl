bitmap(grid::Grid2D) = Bitmap2D(grid)
bitmap(grid::Grid3D) = Bitmap3D(grid)

empty(M::Bitmap2D) = Bitmap2D(M.grid)
empty(M::Bitmap3D) = Bitmap3D(M.grid)
empty(grid::Grid2D) = Bitmap2D(grid)
empty(grid::Grid3D) = Bitmap3D(grid)

grid(B::AbstractBitmap) = B.grid
