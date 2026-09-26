#pragma once
#include "../typedefs/header.hpp"
#include <filesystem>
#include <format>
#include <stdexcept>
#include <string_view>

namespace MathEngine::IO
{ // namespace MathEngine::IO

enum class FPFormat { Fixed, Scientific, Default};
enum class Alignment {Left, Center, Right, None};

struct WriteOptions
{
    std::filesystem::path path;
    std::string_view      separator;
    std::string_view      header;
    std::string_view      comment;
    std::string_view      footer;
    size_t                colWidth;
    int                   precision;
    FPFormat              format;
    Alignment             alignment;
    bool                  append;
    bool                  binary;
};

namespace Details
{
    template <typename T>
    inline std::string BuildFormatSpecifier(const WriteOptions& wOpts)
    {
        std::string fmt = "{:";
        if (wOpts.colWidth>0)
        {
            switch(wOpts.alignment)
            {
                case Alignment::Left:   fmt += '<'; break;
                case Alignment::Center: fmt += '^'; break;
                case Alignment::Right:  fmt += '>'; break;
                case Alignment::None:               break;
                default:                            break;
            }
            fmt += std::to_string(wOpts.colWidth);
        }

        if constexpr (FPNumber<T>)
        {
            if (wOpts.precision>=0) fmt += '.' + std::to_string(wOpts.precision);

            switch(wOpts.format)
            {
                case FPFormat::Fixed:      fmt += 'f'; break;
                case FPFormat::Scientific: fmt += 'e'; break;
                case FPFormat::Default:    fmt += 'e'; break;
                default:                               break;
            }
        }
            
        fmt += '}'; return fmt;
    }

    template <typename T>
    inline void WriteRawBytes(std::ofstream& f, std::span<const T> data)
    {
        std::span<const std::byte> bytes = std::as_bytes(data);

        f.write(reinterpret_cast<const char*>(bytes.data()),bytes.size_bytes());
    }

    inline void WriteHeaderAndComment(std::ofstream& f, const WriteOptions& wOpts)
    {
        if (!wOpts.header.empty())
        {
            f << "# " << wOpts.header << '\n';
            if (!wOpts.comment.empty()) f << "## " << wOpts.comment << '\n';
        }
    }
} // End MathEngine::IO::Details namespace

template <Number T>
inline void WriteMatrix(const Matrix<T>& mat, const WriteOptions& wOpts)
{
    if (wOpts.path.has_parent_path()) std::filesystem::create_directories(wOpts.path.parent_path());

    auto mode = std::ios::out | (wOpts.append?std::ios::app:std::ios::trunc);
    if (wOpts.binary) mode |= std::ios::binary;

    std::ofstream f(wOpts.path,mode);
    if (!f.is_open()) throw std::runtime_error("Could not open file at: "+wOpts.path.string());

    const size_t rows = mat.Rows();
    const size_t cols = mat.Cols();
    if (wOpts.binary)
    {
        Details::WriteRawBytes(f, std::span{&rows,1});
        Details::WriteRawBytes(f, std::span{&cols,1});

        Details::WriteRawBytes(f, std::span{mat.ptr(),rows*cols}); return;
    }

    Details::WriteHeaderAndComment(f, wOpts);

    // const std::string fmt = Details::BuildFormatSpecifier<std::remove_pointer_t<decltype(mat.ptr())>>(wOpts);
    const std::string fmt = Details::BuildFormatSpecifier<T>(wOpts);

    std::string buffer;
    buffer.reserve(rows*cols*(wOpts.colWidth>0?wOpts.colWidth+wOpts.separator.size():20)+(!wOpts.footer.empty()?wOpts.footer.size():100));
    for (size_t r{0}; r<rows; ++r)
    {
        for (size_t c{0}; c<cols; ++c)
        {
            std::vformat_to(std::back_inserter(buffer),fmt,std::make_format_args(mat[r,c]));
            const bool isLast = (c==cols-1);
            if (!isLast) buffer.append(wOpts.separator); 
        }
        buffer.append("\n");
    }

    if (!wOpts.footer.empty())
    {
        buffer.append(wOpts.footer);
        buffer.append("\n");
    }
    f.write(buffer.data(),buffer.size());
}

template <typename T>
inline void WriteVector(std::span<const T> vecView, const WriteOptions& wOpts)
{
    if (wOpts.path.has_parent_path()) std::filesystem::create_directories(wOpts.path.parent_path());

    auto mode = std::ios::out | (wOpts.append?std::ios::app:std::ios::trunc);
    if (wOpts.binary) mode |= std::ios::binary;

    std::ofstream f(wOpts.path,mode);
    if (!f.is_open()) throw std::runtime_error("Could not open file at: "+wOpts.path.string());

    const size_t length = vecView.size();
    if (wOpts.binary)
    {
        Details::WriteRawBytes(f, std::span{&length,1});

        Details::WriteRawBytes(f, vecView); return;
    }
    
    Details::WriteHeaderAndComment(f, wOpts);

    // const std::string fmt = Details::BuildFormatSpecifier<std::remove_pointer_t<decltype(mat.ptr())>>(wOpts);
    const std::string fmt = Details::BuildFormatSpecifier<T>(wOpts);

    std::string buffer;
    buffer.reserve(length*(wOpts.colWidth>0?wOpts.colWidth+wOpts.separator.size():20)+(!wOpts.footer.empty()?wOpts.footer.size():100));
    for (size_t i{0}; i<length; ++i)
    {
        std::vformat_to(std::back_inserter(buffer),fmt,std::make_format_args(vecView[i]));
        const bool isLast = (i==length-1);
        if (!isLast) buffer.append(wOpts.separator); 
    }
    buffer.append("\n");

    if (!wOpts.footer.empty())
    {
        buffer.append(wOpts.footer);
        buffer.append("\n");
    }
    f.write(buffer.data(),buffer.size());
}

} // End MathEngine::IO namespace

