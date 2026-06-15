#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

#include "polypartition.h"

static bool read_polygon_txt(const std::string& path, TPPLPoly& poly) {
    std::ifstream fin(path);
    if (!fin) {
        std::cerr << "ERROR: cannot open input file: " << path << "\n";
        return false;
    }

    long n;
    fin >> n;
    if (!fin || n < 3) {
        std::cerr << "ERROR: input must start with vertex count n >= 3.\n";
        return false;
    }

    poly.Init(n);
    poly.SetHole(false);

    for (long i = 0; i < n; ++i) {
        double x, y;
        fin >> x >> y;
        if (!fin) {
            std::cerr << "ERROR: failed reading vertex " << i << ".\n";
            return false;
        }
        poly[i].x = x;
        poly[i].y = y;
        poly[i].id = static_cast<int>(i);
    }

    poly.SetOrientation(TPPL_ORIENTATION_CCW);
    return true;
}

static bool write_pieces_json(
    const std::string& path,
    const std::string& mode,
    int ok,
    const TPPLPolyList& pieces
) {
    std::ofstream fout(path);
    if (!fout) {
        std::cerr << "ERROR: cannot open output file: " << path << "\n";
        return false;
    }

    fout << std::setprecision(17);
    fout << "{\n";
    fout << "  \"ok\": " << (ok ? "true" : "false") << ",\n";
    fout << "  \"mode\": \"" << mode << "\",\n";
    fout << "  \"piece_count\": " << pieces.size() << ",\n";
    fout << "  \"pieces\": [\n";

    long piece_idx = 0;
    for (TPPLPolyList::const_iterator it = pieces.begin(); it != pieces.end(); ++it, ++piece_idx) {
        const TPPLPoly& p = *it;
        fout << "    [";
        for (long i = 0; i < p.GetNumPoints(); ++i) {
            if (i > 0) fout << ", ";
            fout << "[" << p[i].x << ", " << p[i].y << "]";
        }
        fout << "]";
        if (piece_idx + 1 < static_cast<long>(pieces.size())) fout << ",";
        fout << "\n";
    }

    fout << "  ]\n";
    fout << "}\n";
    return true;
}

int main(int argc, char** argv) {
    if (argc != 4) {
        std::cerr
            << "Usage:\n"
            << "  polypp_cli.exe triangulate_mono input.txt output.json\n"
            << "  polypp_cli.exe convex_hm       input.txt output.json\n";
        return 2;
    }

    std::string mode = argv[1];
    std::string input_path = argv[2];
    std::string output_path = argv[3];

    TPPLPoly poly;
    if (!read_polygon_txt(input_path, poly)) {
        return 3;
    }

    TPPLPartition partition;
    TPPLPolyList pieces;
    int ok = 0;

    if (mode == "triangulate_mono") {
        ok = partition.Triangulate_MONO(&poly, &pieces);
    } else if (mode == "convex_hm") {
        ok = partition.ConvexPartition_HM(&poly, &pieces);
    } else {
        std::cerr << "ERROR: unknown mode: " << mode << "\n";
        return 4;
    }

    if (!write_pieces_json(output_path, mode, ok, pieces)) {
        return 5;
    }

    if (!ok) {
        std::cerr << "ERROR: PolyPartition returned failure for mode: " << mode << "\n";
        return 6;
    }

    std::cout << "OK: " << mode << " produced " << pieces.size() << " pieces.\n";
    return 0;
}

