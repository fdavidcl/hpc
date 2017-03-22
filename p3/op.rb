#!/usr/bin/env ruby
# coding: utf-8

def selections amount
  Enumerator.new do |yielder|
    selected = Array.new(amount, 0)
    yielder << selected
    
    while selected.any? &:zero?
      i = selected.length - 1
      
      while selected[i] != 0
        selected[i] = 0
        i -= 1
      end
      
      selected[i] = 1
      yielder << selected
    end
  end
end

# tabla 0-1 para una operación dada
def table n, &block
  args = selections n
  results = []
  
  begin
    loop do
      aa = args.next
      results << [*aa, block.call(*aa)]
    end
  rescue StopIteration # paramos cuando termine el enumerador
  rescue StandardError => e
    puts e
  end

  results
end

def show_table *args, &block
  puts (table(*args, &block).map { |r| r.join " " }.join("\n"))
end


triangle = ->(left, right, up, down) {
  down * (1 + [left, right, up].min)
}

# size[i][j] = min(size[i - 1][j - 1], min(size[i - 1][j], size[i][j - 1])) + 1;

leetcode_sol = ->(left, right, up, down) {
  if down == 1
    [up, right, left].min + 1
  else
    down
  end
}

show_table 4, &triangle
show_table 4, &leetcode_sol
