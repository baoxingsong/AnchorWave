#pragma once

#include <string>

namespace anchorwave {
namespace io {

// Inspect the file signature rather than its suffix. BGZF uses the gzip magic
// bytes and is therefore detected here as well.
bool isGzipCompressed(const std::string &path);

// Return path unchanged for plain input. Gzip/BGZF input is decompressed once
// into a process-lifetime temporary file so existing seek/pread based FASTA
// indexing remains valid. The temporary file is removed on normal process exit.
std::string materializeInputFile(const std::string &path);

}  // namespace io
}  // namespace anchorwave
