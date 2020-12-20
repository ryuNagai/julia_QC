function test()
    a = 0
    for i in 1:10000
        a = a + 1
    end
    println(a)
end

@time test()