package com.oktsrl.utils;

import java.io.Serializable;

public final class Action implements Comparable<Action>,Serializable{

	private static final long serialVersionUID = -1013389969348628100L;

	public int user;
    public int item;
    public double rating;
    public long timestamp;
    
    
    public long getTimestamp() {
        return timestamp;
    }

    public void setTimestamp(long timestamp) {
        this.timestamp = timestamp;
    }

    public int getUser() {
        return user;
    }

    public int getItem() {
        return item;
    }

    public double getRating() {
        return rating;
    }

    @Override
    public String toString() {
        return String.format("action: user=%d, item=%d, rating=%.1f", user, item , rating);
    }

    public Action() {
    }

    public Action(int user, int item, double rating) {
        this.user = user;
        this.item = item;
        this.rating = rating;
        this.timestamp=0L;
    }
    
    public Action(int user, int item, double rating,long timestamp) {
        this.user = user;
        this.item = item;
        this.rating = rating;
        this.timestamp=timestamp;
    }

    
    
    @Override
    public int hashCode() {
        final int prime = 31;
        int result = 1;
        result = prime * result + item;
        result = prime * result + user;
        return result;
    }

    @Override
    public boolean equals(Object obj) {
        if (this == obj)
            return true;
        if (obj == null)
            return false;
        if (getClass() != obj.getClass())
            return false;
        Action other = (Action) obj;
        if (item != other.item)
            return false;
        if (user != other.user)
            return false;
        return true;
    }

    @Override
    public int compareTo(Action o) {
        if(timestamp<o.timestamp)
            return -1;
        if(timestamp>o.timestamp)
            return 1;
        return 0;
    }
    
}