const POSIX_FADV_NORMAL = 0      # No further special treatment.
const POSIX_FADV_RANDOM = 1      # Expect random page references.
const POSIX_FADV_SEQUENTIAL =  2 # Expect sequential page references.
const POSIX_FADV_WILLNEED = 3    # Will need these pages.
const POSIX_FADV_DONTNEED = 4    # Don't need these pages.
const POSIX_FADV_NOREUSE  = 5    # Data will be accessed once.

if is_linux()
  function fadvise(fd, offset, len, advice)
    ccall(:posix_fadvise, Cint, (Cint, Cint, Cint, Cint), fd, offset, len, advice)
  end
else
  function fadvise(fd, offset, len, advice)
    warn("Cache eviction is not yet supported on this platform")
  end
end

function uncache(filename :: String)
  file = open(filename, "r")
  uncache(file)
  close(file)
end

function uncache(io :: IOStream)
  fadvise(fd(io), 0, 0, POSIX_FADV_DONTNEED)
end

