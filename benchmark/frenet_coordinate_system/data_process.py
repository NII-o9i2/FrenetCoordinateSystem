


if __name__ == "__main__":
    # data vx 
    data_x_file = open("data_x.txt")
    match_string = "frenet_data vx"
    vx = []
    for line in data_x_file.readlines():
        for cnt in range(len(line)):
            str_tmp = ""
            max_len = int(min(cnt + 14,len(line)-1))
            if line[int(cnt):max_len] == match_string:
                str_tmp = line[max_len:len(line)-1]
                vx.append(float(str_tmp))
                break
    print("len vx" , len(vx))
    print(vx)
    data_x_file.close()
    
    # data vy
    data_y_file = open("data_y.txt")
    match_string = "frenet_data vy"
    vy = []
    for line in data_y_file.readlines():
        for cnt in range(len(line)):
            str_tmp = ""
            max_len = int(min(cnt + 14,len(line)-1))
            if line[int(cnt):max_len] == match_string:
                str_tmp = line[max_len:len(line)-1]
                vy.append(float(str_tmp))
                break
    print("len vy" , len(vy))            
    print(vy)
    data_x_file.close() 