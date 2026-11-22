class Cellcoal < Formula
  desc "Coalescent simulator of single-cell NGS data (somatic mutations, ADO, doublets, reads)"
  homepage "https://github.com/dapogon/cellcoal"
  url "https://github.com/dapogon/cellcoal/archive/refs/tags/v2.6.1.tar.gz"
  sha256 "bda773b5fbc2968cec6222fd56b600f53d806299d9d046795c0a914d94f99c18"
  license "GPL-2.0-or-later"

  depends_on "libomp"

  def install
    system "make"
    system "make", "install", "PREFIX=#{prefix}"
  end

  test do
    system bin/"cellcoal", "-n1", "-s10", "-l200", "-u1e-6", "-e10000", "-b1", "-#1"
    assert_predicate testpath/"Results", :directory?
  end
end
