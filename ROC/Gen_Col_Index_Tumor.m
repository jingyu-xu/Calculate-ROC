function col_vec_tumor = Gen_Col_Index_Tumor(z,x,y,s1_ss,s2_ss,diameter,k)

col_vec_tumor = zeros(1,diameter^3);
for zi = z:z+diameter-1
    for xi = x:x+diameter-1
        for yi = y:y+diameter-1
            col_vec_tumor(1,k) = (zi-1)*(s1_ss*s2_ss)+(xi-1)*s2_ss+yi;
            k = k+1;
        end
    end
end
        