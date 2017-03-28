#!/usr/bin/env ruby
# coding: utf-8

if ARGV.length < 2
  raise ArgumentError, "Necesito un archivo de matriz y un número de particiones como parámetro"
end

filename = ARGV[0].strip
partitions = ARGV[1].strip.to_i

lines = File.readlines(filename)

slice_length = lines.length/partitions
lines.each_slice(slice_length).each_with_index do |partition, i|
  File.write("#{filename}-#{partitions}p#{i}", partition.join)
end
