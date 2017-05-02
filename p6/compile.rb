#!/usr/bin/env ruby

CU = "nvcc"
SRC = "toy.cu io.c"
BASENAME = "toy_"

CFLAGS   = "-march=native -Ofast -Wall -Wl,--no-as-needed"
CUFLAGS  = "-D_MWAITXINTRIN_H_INCLUDED -I./"

def compile exe, params
  paramline = params.map do |key, value|
    value ? "-D#{key}=#{value}" : ""
  end.join " "

  compilation = "#{CU} #{CUFLAGS} -o #{exe} #{SRC} #{paramline}"
  puts compilation
  system compilation
end

DEV_ID = 0

def execute exe, num_ops, input, output
  `./#{exe} #{DEV_ID} #{num_ops} #{input} #{output}`
end

PARAMS = {
  TYPE: ["float", "int", "double"],
  OP: [false, "SQUARE", "DIVISION"],
  UNSOLOSTREAM: [false],
}

combinations = PARAMS[:TYPE].product(PARAMS[:OP]).product(PARAMS[:UNSOLOSTREAM])

EXE_PARAMS = {
  num_op: [1, 10, 100, 1000],
  input: Dir["data/*"]
}

data = combinations.map do |pars|
  pars.flatten!
  exe = "toy_#{pars.join "_"}"
  compile exe, PARAMS.keys.zip(pars).to_h
  exe

  EXE_PARAMS[:num_op].product(EXE_PARAMS[:input]).map do |num_op, input|
    output_m = "out/#{exe}_#{input}_#{num_op}.rua"
    out = execute exe, num_op, input, output_m
    "#{pars.join ", "}, #{num_op}, #{input}, #{out}"
  end
end.flatten

header = "type, op, unstream, num_op, matriz, t_acceso, creacion, ejecucion, computacion, glopbyte, gflops\n"

File.write("p6.csv", header + data.join("\n"))
