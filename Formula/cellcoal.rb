class Cellcoal < Formula
  desc "Coalescent simulator of single-cell NGS data (somatic mutations, ADO, doublets, reads)"
  homepage "https://github.com/dapogon/cellcoal"
  url "https://github.com/dapogon/cellcoal/archive/refs/tags/v2.7.4.tar.gz"
  sha256 "259c95e77182dd1751717eeafd01a3d8eda6a8eb7031b58b71424bbd5e5861b2"
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
